function d = compute_scores(X, dtype, dparams)
%COMPUTE_SCORES Computes a score/probability for some data
%   [d] = COMPUTE_SCORES(X, dtype, dparams)
%   computes some score/probability of X based on the given parameters.
%
%   X should be a matrix of size N x D, having N D-dimensional observations
%   for which the scores will be computed.
%
%   dtype specifies the type of discriminant to be used. Currently defined
%   are:
%   - 'dlinear'
%   - 'dlinearnocov'
%   - 'dlog'
%   - 'dmvnpdf'
%   - 'dmvnpdfnorm'
%   - 'dmvnpdfnormspool'
%   - 'dmdpdf'
%
%   dparams is a structure that provides the necessary parameters for
%   computing the discriminant, including the descriptive statistics (mean,
%   covariance). These paramters include:
%   - 'values': a vector of length nr_values specifying the values
%   corresponding to the mean vectors provided in this structure.
%   - 'mvec': a matrix of size nr_values x nr_samples, with the mean
%   vectors of each group (value). This parameter is mandatory.
%   - 'scov': covariance matrix, given for each value or a single (pooled)
%   matrix. Size is nr_samples x nr_samples.
%   - 'sinv': the inverse of the covariance matrix, given for each value or
%   a single matrix (pooled covariance). Size is nr_samples x nr_samples.
%   - 'slogdet': log-determinant of covariance matrix for each value/group.
%   - 'plaintext': used only with dtype='dpa', to provide the corresponding
%   plaintext of each attack trace. This should be a matrix of size
%   nr_traces x nr_groups. Currently not implemented.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)
%
%   See the "Applied Multivariate Statistical Analysis" book for details.
%
%   See also evaluate_discriminant, evaluate_mdlike, mvnlikelog, mdlike.


%% Check and initialise parameters
if ~isfield(dparams, 'mvec')
    error('Missing mvec parameter');
end
N = size(X, 1);
D = size(X, 2);
nr_groups = size(dparams.mvec, 1);
if size(dparams.mvec, 2) ~= D
    error('Incompatible size of mvec with X');
end
d = zeros(nr_groups, 1);

%% Compute scores from given parameters
if strcmp(dtype, 'dlinear')
    %% Compute linear discriminant using common covariance
    ct = -N/2;
    xs = ones(1,N) * X;
    for k = 1:nr_groups
        d(k) = dparams.mvec(k,:)*dparams.sinv*xs' + ...
               ct*(dparams.mvec(k,:)*dparams.sinv*dparams.mvec(k,:)');        
    end
elseif strcmp(dtype, 'dlinearnocov')
    %% Compute linear discriminant using identity covariance
    ct = -N/2;
    xs = ones(1,N) * X;
    for k = 1:nr_groups
        d(k) = dparams.mvec(k,:)*xs' + ...
               ct*(dparams.mvec(k,:)*dparams.mvec(k,:)');        
    end
elseif strcmp(dtype, 'dlog')
    %% Compute log-likelihood discriminant
    ct1 = -N/2;
    ct2 = -1/2;
    for k = 1:nr_groups
        dsum = 0;
        for j = 1:N
            x = X(j,:) - dparams.mvec(k,:);
            dsum = dsum + x*dparams.sinv(:,:,k)*x';
        end
        d(k) = ct1*dparams.slogdet(k) + ct2*dsum;
    end
elseif strcmp(dtype, 'dmvnpdf')
    %% Compute multivariate normal pdf, going first in log domain
    ct1 = -(D*N/2)*log(2*pi);
    ct2 = -N/2;
    ct3 = -1/2;
    for k = 1:nr_groups
        dsum = 0;
        for j = 1:N
            x = X(j,:) - dparams.mvec(k,:);
            dsum = dsum + x*dparams.sinv(:,:,k)*x';
        end
        vlog = ct1 + ct2*dparams.slogdet(k) + ct3*dsum;
        d(k) = exp(vlog);
    end   
elseif strcmp(dtype, 'dmvnpdfnorm') || strcmp(dtype, 'dmvnpdfnormspool')
    %% Compute mvn pdf, normalising all probas after each attack trace
    for k=1:nr_groups
        d(k) = 1;
    end
    
    for j = 1:N
        for k = 1:nr_groups
            pmvn = mvnpdfinv(X(j,:), dparams.mvec(k,:), ...
                dparams.sinv(:,:,k), dparams.slogdet(k));
            d(k) = d(k) * pmvn;            
        end
        sdk = sum(d);        
        for k = 1:nr_groups
            d(k) = d(k) / sdk;
        end
    end    
elseif strcmp(dtype, 'dmdpdf')
    if size(dparams.sinv, 3) > 1
        error('dmdpdf supported for single pooled covariance');
    end
    %% Compute Mahalanobis distance pdf based on Chi-square PDF
    % Compute first in log domain to avoid underflow for large N
    for k = 1:nr_groups
        vlog = 0;
        for j = 1:N
            x = X(j,:) - dparams.mvec(k,:);
            md = x*dparams.sinv*x';
            vlog = vlog + log(chi2pdf(md, D));
        end
        d(k) = exp(vlog);
    end    
end

end



