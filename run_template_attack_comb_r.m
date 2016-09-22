function [results] = run_template_attack_comb_r(...
    X_profile, V_profile, ...
    X_attack, X_attack_input, ...
    V_attack, V_key, comb_func, ...
    cmethod, cparams, ...
    dtype, ...
    nr_traces_vec, rand_iter, eparams)
%RUN_TEMPLATE_ATTACK_COMB Runs a template attack against a crypto function
%   [results] = RUN_TEMPLATE_ATTACK_COMB(...
%       X_profile, V_profile, ...
%       X_attack, X_attack_input, ...
%       V_attack, V_key, comb_func, ...
%       cmethod, cparams, ...
%       dtype, ...
%       nr_traces_vec, eparams)
%   runs a template attack on some given combination function (e.g. AES
%   S-box), using data given as matrices for profile and attack,
%   and returns a results scores/probabilities of each possible value/group.
%
%   X_profile should be a matrix of size nr_trials x nr_samples x nr_values,
%   where:
%   - nr_trials is the number of traces each having nr_samples leakage
%     samples.
%   - nr_values is the number of groups for which to run the profiling and
%   more generally the template attack.
%   This matrix should contain the data for profiling, with data grouped by
%   target values.
%
%   V_profile should be a vector of length nr_values, containing the
%   target value (e.g. S-box(key XOR plaintext)) for each group of leakage
%   traces given in X_profile.
%
%   X_attack should be a matrix of size nr_trials x nr_samples x nr_bytes
%   containing the leakage traces for attacking different target values
%   (one nr_trials x nr_samples matrix per byte value). This method
%   will compute the discriminants/scores per byte value.
%
%   X_attack_input should be a matrix of size nr_trials x nr_bytes
%   specifying the plaintext used for each of the leakage traces in the
%   X_attack matrix. During the template attack, this method will compute
%   the value of the combination function (see below) for each key byte
%   hypopthesis.
%
%   V_attack should be a vector of length nr_bytes specifying the actual
%   values of the target bytes for which the matrix X_attack was given.
%
%   V_key should be a vector of length nr_keys specifying all the possible
%   values of the key, i.e. this should be the vector of key hypothesis,
%   that will be used during the attack to find the most likely value.
%
%   comb_func should be a method, taking two parameters: a vector of
%   plaintext values and a vector of key bytes, and which should return the
%   target combination of these (e.g. Sbox(p XOR k)). This method should be
%   the same that was used to compute the values in V_profile, and
%   will be used during the attack to compute the output when trying all
%   possible key byte hypothesis. This combination function should expect
%   values starting with 0.
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
%   dtype specifies the type of discriminant to be used. See compute_scores
%   for available types.
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used for each element.
%
%   rand_iter should be a positive integer specifying the number of
%   iterations to run the evaluation (with different attack traces).
%
%   eparams is a structure of extra parameters that may be needed for some
%   options. Some of these extra parameters include:
%   - 'save_xdata': if this field exists and is non zero then the x_profile
%   and x_attack data will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'save_eval': if this field exists and is non zero then the
%   parameters for evaluation (dparams) will be saved. These take
%   considerable amount of space.
%   - 'save_ssp': use this to save the SSP matrices M,B,W.
%   - 'v_attack': if this field is given then it should contain a vector of
%   the same length this method's parameter V_attack. This vector will
%   replace V_attack. This may be useful to provide V_attack as values
%   covering a wide range (e.g. 16-bit data) but then running the
%   evaluation on only a part of the data (e.g. the most significant 8
%   bits). I used it to test attacks on 8-bit data influenced by pipeline
%   of other 8-bit data. Hence the original V_attack had 16-bit data but
%   then I run the attack on only one byte, assuming the remaining byte is
%   noise that the attack has to deal with.
%   - 'cf_params': parameters for the combination function.
%   - 'use_fa': use factor analysis for covariance matrix.
%   - 'nr_factors': number of factors to use with factor analysis.
%
%   The 'results' structure contains the following:
%   -> results.dtype: the discriminant type string.
%   -> results.dparams: parameters for evaluation (optional).
%   -> results.nr_traces_vec: the vector with number of attack traces.
%   -> results.disc_info: discriminant scores/probabilities, as returned by
%   get_disc_info.
%   -> results.error: an optional messsage in case of an error. This may
%   avoid the abortion of a larger set of runs while allowing to detect the
%   cases in which an error occurred.
%
%   See the paper "Efficient Template Attacks", Choudary and Kuhn, CARDIS 2013.
%
%   See also get_mmap.

%% Check and initialise parameters
n_profile = size(X_profile, 1);
nr_values = length(V_profile);
if nr_values ~= size(X_profile, 3)
    error('Incompatible number of values for X_profile and V_profile');
end
results = [];
results.dtype = dtype;
results.nr_traces_vec = nr_traces_vec;
results.rand_iter = rand_iter;
results.v_profile = V_profile;
results.v_attack = V_attack;
results.v_key = V_key;
results.comb_func = comb_func;
results.cparams = cparams;
results.cmethod = cmethod;
if nargin < 9
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_template_attack_comb_r() ...\n');

%% Get the sums of squares and cross products
fprintf('Obtaining sums of squares and cross products...\n');
[M, B, W] =  compute_ssp(X_profile);
if isfield(eparams, 'save_ssp') && (eparams.save_ssp ~= 0)
    results.M = M;
    results.B = B;
    results.W = W;
end
xmm = mean(M, 1);
toc

%% Estimate correlation via factor analysis if specified
if isfield(eparams, 'use_fa')
    if ~isfield(eparams, 'nr_factors')
        error('Need eparams.nr_factors for famvn');
    else
        nr_factors = eparams.nr_factors;
    end
    C = W / (nr_values*(n_profile - 1));
    [U, S, ~] = svd(C);
    d = diag(S);
    L = U(:,1:nr_factors)*diag(sqrt(d(1:nr_factors)));
    results.L = L;
    P = C - L*L';
    P = diag(P);
    results.P = P;
    CE = L*L' + diag(P); % Estimated covariance matrix from factor analysis        
    W = CE * (nr_values*(n_profile - 1)); % Replace W with data from CE
end

%% Get compression parameters
if strcmp(cmethod, 'sample')
    fprintf('Computing selection curves and selected samples ...\n');
    [curves] = get_signal_strength_ssp(M, B, W, n_profile);
    
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
    Spool = W / (nr_values*(n_profile-1));
    [A ,~, K] = compute_params_lda(B, Spool, nr_values, cparams.lda_threshold);
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

%% Compress/process profile data as necessary
fprintf('Processing profiling data...\n');
X_profile_c = compute_features_data(X_profile, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_profile = X_profile_c;
end
toc

%% Compress/process attack data as necessary
fprintf('Processing attack data...\n');
X_attack_c = compute_features_data(X_attack, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5);
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_attack = X_attack_c;
end
toc

%% Compute templates and evaluation parameters
fprintf('Computing template and evaluation parameters...\n');
[tmiu, tsigma] = compute_template(X_profile_c);
dparams = [];
dparams.values = V_profile;
dparams.mvec = tmiu;
dparams.scov = tsigma;
dparams.plaintext = X_attack_input;
dparams.comb_func = comb_func;
if isfield(eparams, 'cf_params')
    dparams.cf_params = eparams.cf_params;
end

% TODO: implement/add for linearfast
if (strcmp(dtype, 'dlinear') || strcmp(dtype, 'dmdpdf'))
    if isfield(eparams, 'use_fa')
        if strcmp(cmethod, 'LDA') || strcmp(cmethod, 'PCA')
            c0 = pp1' * CE * pp1;
        elseif strcmp(cmethod, 'sample')
            c0 = CE(pp1,pp1);
        end
    else
        c0 = mean(tsigma, 3);
    end
    ic0 = inv(c0);
    dparams.sinv = ic0;
elseif (strcmp(dtype, 'dlog') || strcmp(dtype, 'dmvnpdf'))
    n = size(tsigma, 3);
    tsinv = zeros(size(tsigma));
    tlogdet = zeros(n, 1);   
    for k=1:n
        if isfield(eparams, 'use_fa')
            CK = cov(X_profile(:,:,k));
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
        end
        tsinv(:,:,k) = inv(tsigma(:,:,k));
        tlogdet(k) = logdet(tsigma(:,:,k), 'chol');
    end
    dparams.sinv = tsinv;
    dparams.slogdet = tlogdet;
end

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.dparams = dparams;
end

%% Compute the discriminant scores/probabilities
fprintf('Computing discriminant info...\n');
[results.disc_info] = get_disc_info_comb_r(...
                               X_attack_c, V_attack, V_key, ...
                               nr_traces_vec, rand_iter, ...
                               dtype, dparams, comb_func);
toc

end
