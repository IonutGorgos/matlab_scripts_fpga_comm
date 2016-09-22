function [results] = run_template_attack(...
    m_data_profile, metadata_profile, idx_profile, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams)
%RUN_TEMPLATE_ATTACK Runs a template attack
%   [results] = RUN_TEMPLATE_ATTACK(...
%       m_data_profile, metadata_profile, idx_profile, ...
%       m_data_attack, metadata_attack, idx_attack, ...
%       bytes, atype, cmethod, cparams, discriminant, ...
%       rand_iter, nr_traces_vec, eparams)
%   runs a template attack with the given parameters and returns a results
%   structure that is defined below.
%
%   This method superseeds run_template_attack_e2, as it takes a general
%   memory mapped object which should represent any kind of leakage, as
%   long as it contains a structure similar with the one used in the E2
%   experiments. See also get_mmap. In addition this method should be much
%   faster when either using the same data for profiling and attack or
%   using the same data for multiple consecutive calls of this method.
%
%   m_data_profile and metadata_profile should be the memory mapped object
%   and associated metadata info for the profiling data. Use get_mmap 
%   on the selected data to obtain these objects.
%
%   idx_profile should be a vector of indices specifying which traces from
%   each group should be used for the profile data.
%
%   m_data_attack and metadata_attack should be the memory mapped object
%   and associated metadata info for the attack data. Use get_mmap 
%   on the selected data to obtain these objects. You can pass the same
%   objects for profiling and attack, which should make the attack faster.
%
%   idx_attack should be a vector of indices specifying which traces from
%   each group should be used for the attack data.
%   
%   bytes should be a vector of indices specifying which bytes (starting
%   from 0) will be used for the attack. This might be useful in order to
%   restrict the attack only to the bytes 0-15 for example (i.e. using 4
%   bits).
%
%   atype should be a string specifying the type of template attack to be
%   used. Currently supported are:
%   - 'mvn': which relies on the multivariate normal probability density
%   function to compute templates and probabilities.
%   - 'famvn': similar to 'mvn' but uses factor analysis to estimate the
%   correlation between samples.
%
%   cmethod should be a string specifying the compression method. Currently
%   supported methods are: 'sample', 'PCA' and 'LDA'.
%
%   cparams should be a structure of params specific to the compression
%   method. For each compression method the params are as follows:
%   - 'sample':
%       -> cparams.curve is a string specifying the signal strength curve to
%       be used.
%       -> cparams.sel is a string specifying the class of selection.
%       -> cparams.p1 is a parameter for the class of selection.
%       -> cparams.p2 is an additional optional parameter.
%   - 'PCA':
%       -> cparams.pca_threshold
%       -> cparams.pca_alternate
%       -> cparams.pca_dimensions
%   - 'LDA':
%       -> cparams.lda_threshold
%       -> cparams.lda_dimensions
%
%   discriminant should be a string specifying the type of discriminant to
%   be used. The possible options are:
%   - 'linear': uses a pooled common covariance matrix with a linear
%      discriminant.
%   - 'linearnocov': does not use a covariance matrix. Might be useful in
%      particular with LDA, where the covariance should be the
%      identity if the eigenvectors are chosen carefully.
%   - 'log': uses individual covariances and log-determinants to compute
%      the group specific templates (mean and covariance).
%
%   rand_iter should be a positive integer specifying the number of
%   iterations to run the evaluation (guessing_entropy) computation. The
%   returned results may contain either the individual or the average
%   results. Check below for details.
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used for each element.
%
%   eparams is a structure of extra parameters that may be needed for some
%   options. Some of these extra parameters include:
%   - 'nr_factors': number of factors to use for factor analysis (which is
%   enabled using the 'famvn' atype).
%   - 'save_xdata': if this field exists and is non zero then the x_profile
%   and x_attack data will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'save_eval': if this field exists and is non zero then the
%   handle_eval pointer and the data for evaluation (generally the
%   templates) will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'use_eig': for LDA, use eigenvectors rather than singular values
%   (which is the preferred method for both PCA and LDA).
%
%   The 'results' structure contains the following:
%   -> results.metadata_profile: the metadata structure for profile.
%   -> results.idx_profile: the idx_profile vector.
%   -> results.metadata_attack: the metadata structure for attack.
%   -> results.idx_attack: the idx_attack vector.
%   -> results.bytes: the bytes vector.
%   -> results.atype: the atype string.
%   -> results.cmethod: the compression method string.
%   -> results.cparams: the cparams structure.
%   -> results.discriminant: the discriminant string.
%   -> results.rand_iter: the number of iterations.
%   -> results.nr_traces_vec: the vector with number of attack traces.
%   -> results.M: the matrix of group means.
%   -> results.B: the between-groups matrix.
%   -> results.W: the matrix of variances and covariances across all data.
%   -> results.L: matrix of factor loadings (only for 'famvn')
%   -> results.P: matrix of individual factors (only for 'famvn')
%   -> results.x_profile: the profiling data, after compression. (optional)
%   -> results.x_attack: the attack data, after compression. (optional)
%   -> results.handle_prepare: function used to extract features for templates.
%   -> results.pp1 ... results.pp5: parameters of results.handle_prepare.
%   -> results.handle_eval: function used to evaluate templates. (optional)
%   -> results.pe3 ... results.pe6: parameters for results.handle_eval. (optional)
%   -> results.success_info: guessing entropy information, as returned by
%   the get_success_info_like method.
%
%   See the paper "Efficient Template Attacks" by Omar Choudary and Markus
%   Kuhn, presented at CARDIS 2013.

%% Check and initialise parameters
np = length(idx_profile);
nr_groups = length(bytes);
results.metadata_profile = metadata_profile;
results.idx_profile = idx_profile;
results.metadata_attack = metadata_attack;
results.idx_attack = idx_attack;
results.bytes = bytes;
results.atype = atype;
results.cmethod = cmethod;
results.cparams = cparams;
results.discriminant = discriminant;
results.rand_iter = rand_iter;
results.nr_traces_vec = nr_traces_vec;
if nargin < 14
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_template_attack() ...\n');

%% Get the sums of squares and cross products
fprintf('Obtaining sums of squares and cross products...\n');
[M, B, W] = compute_ssp_e2_mmap(m_data_profile, metadata_profile, ...
                                idx_profile, bytes);
results.M = M;
results.B = B;
results.W = W;
xmm = mean(M, 1);
toc

%% Estimate correlation via factor analysis if 'famvn' specified
if strcmp(atype, 'famvn')
    if ~isfield(eparams, 'nr_factors')
        error('Need eparams.nr_factors for famvn');
    else
        nr_factors = eparams.nr_factors;
    end
    C = W / (nr_groups*(np - 1));
    [U, S, ~] = svd(C);
    d = diag(S);
    L = U(:,1:nr_factors)*diag(sqrt(d(1:nr_factors)));
    results.L = L;
    P = C - L*L';
    P = diag(P);
    results.P = P;
    CE = L*L' + diag(P); % Estimated covariance matrix from factor analysis        
    W = CE * (nr_groups*(np - 1)); % Replace W with data from CE
end

%% Get compression parameters
if strcmp(cmethod, 'sample')
    fprintf('Computing selection curves and selected samples ...\n');
    [curves] = get_signal_strength_ssp(M, B, W, np);
    
    if ~isfield(cparams, 'p2')
        cparams.p2 = [];
    end
    interest_points = get_selection(curves.(cparams.curve), cparams.sel, ...
        cparams.p1, cparams.p2);
    
    handle_prepare = @prepare_data_template_xmm;
    pp1 = interest_points;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(cmethod, 'PCA')
    fprintf('Computing PCA parameters...\n');
    [U, ~, xmm, K] = compute_params_pca(M, cparams.pca_threshold, ...
                                        cparams.pca_alternate);
    if isfield(cparams, 'pca_dimensions') && (cparams.pca_dimensions > 0)
        U = U(:,1:cparams.pca_dimensions);
    else
        U = U(:,1:K);
    end
    
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = U;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(cmethod, 'LDA')
    fprintf('Computing Fishers LDA parameters...\n');
    Spool = W / (nr_groups*(np-1));
    if ~isempty(eparams) && eparams.use_eig
        use_eig = 1;
    else
        use_eig = 0;
    end
    [A ,~, K] = compute_params_lda(B, Spool, nr_groups, cparams.lda_threshold, use_eig);
    if isfield(cparams, 'lda_dimensions') && (cparams.lda_dimensions > 0)
        FW = A(:,1:cparams.lda_dimensions);        
    else
        FW = A(:,1:K);
    end

    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = FW;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
else
    error('Unknown compression method: %s', cmethod);
end

%% Store handle_prepare data
results.handle_prepare = handle_prepare;
results.pp1 = pp1;
results.pp2 = pp2;
results.pp3 = pp3;
results.pp4 = pp4;
results.pp5 = pp5;

%% Load raw leakage data for profile
% Using evalc to avoid console output
fprintf('Computing profiling data...\n');
[~, x_profile] = evalc(['compute_features_e2_mmap(', ...
                        'm_data_profile, metadata_profile, idx_profile,', ...
                        'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_profile = x_profile;
end
toc

%% Load raw leakage data for attack
fprintf('Computing attack data...\n');
[~, x_attack] = evalc(['compute_features_e2_mmap(', ...
                        'm_data_attack, metadata_attack, idx_attack,', ...
                        'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_attack = x_attack;
end
toc

%% Compute templates
if strcmp(atype, 'mvn')
    fprintf('Computing mvn template and evaluation parameters...\n');
    [tmiu, tsigma] = compute_template(x_profile);
    
    handle_eval = @evaluate_discriminant;
    if strcmp(discriminant, 'linear')
        c0 = mean(tsigma, 3);
        ic0 = inv(c0);
        pe3 = tmiu;
        pe4 = ic0;
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'linearnocov')
        pe3 = tmiu;
        pe4 = [];
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'log')
        n = size(tsigma, 3);
        tsinv = zeros(size(tsigma));
        tlogdet = zeros(n, 1);
        for k = 1:n
            tsinv(:,:,k) = inv(tsigma(:,:,k));
            tlogdet(k) = logdet(tsigma(:,:,k), 'chol');
        end
        pe3 = tmiu;
        pe4 = tsinv;
        pe5 = tlogdet;
        pe6 = [];
    else
        error('discriminant not supported for mvn: %s', discriminant);
    end 
elseif strcmp(atype, 'famvn')
    fprintf('Computing mvn template and evaluation parameters...\n');
    [tmiu, tsigma] = compute_template(x_profile);
    
    handle_eval = @evaluate_discriminant;
    if strcmp(discriminant, 'linear')
        % transform full pooled matrix obtained from factor analysis
        if strcmp(cmethod, 'LDA') || strcmp(cmethod, 'PCA')
            c0 = pp1' * CE * pp1;
        elseif strcmp(cmethod, 'sample')
            c0 = CE(pp1,pp1);
        end
        ic0 = inv(c0);
        pe3 = tmiu;
        pe4 = ic0;
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'linearnocov')
        pe3 = tmiu;
        pe4 = [];
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'log')
        % compute each individual covariance matrix, then apply factor
        % analysis and transform as selected by the compression method.
        % This is not optimal, given that I already compute templates in
        % compute_features_e2_mmap, but should work fine for now.
        % However, note that this may take a while, in particular the SVD
        % computation for each of the covariances.
        n = size(tsigma, 3);
        tsinv = zeros(size(tsigma));
        tlogdet = zeros(n, 1);
        for k=1:n
            kindex = find(m_data_profile.data(1).B(2,:)==bytes(k));
            lindex = kindex(idx_profile);
            % Below note transpose and conversion in case we had integer class
            X = double(m_data_profile.data(1).X(:,lindex)');
            CK = cov(X);
            [UK, SK, ~] = svd(CK);
            dk = diag(SK);
            LK = UK(:,1:nr_factors)*diag(sqrt(dk(1:nr_factors)));
            PK = CK - LK*LK';
            PK = diag(PK);
            CKE = LK*LK' + diag(PK); % Estimated covariance matrix from factor analysis
            if strcmp(cmethod, 'LDA') || strcmp(cmethod, 'PCA')
                tsigma(:,:,k) = pp1' * CKE * pp1;
            elseif strcmp(cmethod, 'sample')
                tsigma(:,:,k) = CKE(pp1,pp1);
            end
            tsinv(:,:,k) = inv(tsigma(:,:,k));
            tlogdet(k) = logdet(tsigma(:,:,k), 'chol');
        end
        pe3 = tmiu;
        pe4 = tsinv;
        pe5 = tlogdet;
        pe6 = [];    
    else
        error('discriminant not supported for famvn: %s', discriminant);
    end
else
    error('template attack type not supported: %s', atype);
end

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.handle_eval = handle_eval;
    results.pe3 = pe3;
    results.pe4 = pe4;
    results.pe5 = pe5;
    results.pe6 = pe6;
end

%% Compute the success information
fprintf('Computing success info...\n');
[results.success_info] = get_success_info_like(x_attack, rand_iter, ...
                                       nr_traces_vec, ...
                                       handle_eval, pe3, pe4, pe5, pe6);
toc

end

