function [U, D, xmm, K] = compute_params_pca(M, threshold, alternate)
%COMPUTE_PARAMS_PCA Compute PCA parameters
%   [U, D, xmm, K] = COMPUTE_PARAMS_PCA(M, threshold, alternate)
%   computes the PCA parameters (U, D, xmm, K) from the given data.
%
%   This method uses principal component analysis (PCA) to compute the
%   matrix of eigenvectors U and eigenvalues D from the given data in M.
%
%   M must be a matrix of size nr_groups x nr_trace_points, containing
%   the precomputed mean trace values for each group.
%
%   threshold is an optional parameter specifying the threshold used to
%   select the first K dimensions for dimensionality reduction using the
%   cumulative percentage of total variation. Pass [] (empty) to use the
%   default of 95.
%
%   alternate should be passed True if the alternative method proposed by
%   Standaert et al. should be used. Note that this is only useful if the
%   first dimension is small.
%
%   This method returns:
%   - the matrix U containing the eigenvectors. Use the vector D
%     or the output K to determine how many to use in analysis.
%   - the vector D containing the lambda values that correspond to each
%     eigenvector (D is basically the diagonal of the S matrix returned
%     by the SVD algorithm).
%   - the vector xmm containing the average of all the mean traces. This
%     vector can be used to normalize input data before projection.
%   - an optional number of components K that represent the number of
%     eigenvalues needed to reach the specified threshold of the total
%     variance. This is the cummulative variance method of determining the
%     number of components to use.
%
%   See also the paper: "Template Attacks in Principal Subspaces",
%   by Archambeau et al.

%% Initialize and check parameters
nr_groups = size(M, 1);
if nargin < 2 || isempty(threshold)
    threshold = 0.95;
end
if nargin < 3 || isempty(alternate)
    alternate = false;
end

%% Compute the PCA parameters
xmm = mean(M, 1);
X = M - repmat(xmm, nr_groups, 1);

if alternate
    % Using Standaert's variant for large data
    [UU, S, ~] = svd((1/nr_groups) * (X * X'));
    SIQ = inv(sqrt(S));
    U = (1/sqrt(nr_groups)) * (X' * UU) * SIQ;
else
    % Use directly SVD which works fine even for n=2500 (takes ~ 14 s)
    [U, S, ~] = svd((1/nr_groups) * (X' * X));
end

% Store the eigenvalues in D
D = diag(S);

%% Return K if requested
if nargout > 3
    for k=1:nr_groups
        f = sum(D(1:k)) / sum(D);
        if f >= threshold
            K = k;
            break
        end
    end
end

end