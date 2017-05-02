function [tmiu, tsigma] = compute_template(X, params)
%COMPUTE_TEMPLATE Compute template parameters
%   [tmiu, tsigma] = COMPUTE_TEMPLATE(X, params)
%   computes the templates (mean, covariance matrix) for the given data.
%
%   X should be a 3-dimensional matrix of size
%   nr_samples x nr_interest_points x nr_groups containing the data for
%   which templates should be computed. The templates (tmiu, tsigma) are
%   computed for each group. That is, for group k the templates will be
%   computed from X(:,:,k).
%
%   params is an optional structure containing parameters for this method.
%   The available parameters are:
%   - use_mle: Specify this parameter if the templates
%   should use the maximum likelihood estimate. The maximum
%   likelihood estimate for the covariance matrix divides by nr_samples
%   while the unbiased estimator divides by (nr_samples-1).
%   - use_fa: specify this parameters (e.g. set 'params.use_fa=1') to use
%   factor analysis in the computation of the covariance matrix. If this is
%   used then the following parameter (nr_factors) should also be given.
%   - nr_factors: the number of factors to use for estimation of the
%   covariance.
%
%   This method assumes that any pre-processing has been done and all the
%   given data in X will be used to compute the templates. Use a function
%   such as compute_features to extract only a subset of interest points
%   for each sample trace.
%
%   This method returns the following:
%   - tmiu: a matrix of size nr_groups x nr_interest_points, containing the
%     mean vector of each group.
%   - tsigma: a matrix of size nr_interest_points x nr_interest_points x nr_groups,
%     having the covariance matrix of each group.


%% Initialize and check parameters
nr_samples = size(X,1);
nr_interest_points = size(X, 2);
nr_groups = size(X, 3);
tmiu = zeros(nr_groups, nr_interest_points);
tsigma = zeros(nr_interest_points, nr_interest_points, nr_groups);
mct = 1/nr_samples;
if nargin < 2
    params = [];
end
if ~isfield(params, 'use_mle') 
    sct = (1/(nr_samples-1));
else
    sct = mct;
end
if isfield(params, 'use_fa') && ~isfield(params, 'nr_factors')
    error('For factor analysis please specify number of factors');
end

%% Compute the templates for each group
for k=1:nr_groups
    x = X(:,:,k);
    tmiu(k,:) = mct*(ones(1,nr_samples)*x);
    xm = x - ones(nr_samples,1)*tmiu(k,:);
    tsigma(:,:,k) = sct*(xm'*xm);
    if isfield(params, 'use_fa')
        CK = tsigma(:,:,k);
        [UK, SK, ~] = svd(CK);
        dk = diag(SK);
        LK = UK(:,1:params.nr_factors)*diag(sqrt(dk(1:params.nr_factors)));
        PK = CK - LK*LK';
        PK = diag(PK);
        CKE = LK*LK' + diag(PK);
        tsigma(:,:,k) = CKE;
    end
end

end