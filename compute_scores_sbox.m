function [d, p] = compute_scores_sbox(X, V, dtype, dparams)
%COMPUTE_SCORES_SBOX Computes a score/probability in an AES-SBOX scenario
%   [d, p] = COMPUTE_SCORES_SBOX(X, V, dtype, dparams)
%   computes some score/probability of X based on the given parameters.
%
%   X should be a matrix of size N x D, having N D-dimensional observations
%   for which the scores will be computed.
%
%   V should be a vector of length N (number of traces), where each element
%   is the plaintext corresponding to a trace. For each key hypothesis k*,
%   a value v=aes-sbox(plaintext XOR k*) will be computed and this value
%   will be used for each trace in the computation of the discriminants.
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
%
%   Note 1: that while this method still uses the dlinear computation, there
%   is no complexity advantage in using this over dlog anymore, since we
%   need to use a different mean vector and covariance matrix for each v.
%
%   Note 2: this method assumes that the key hypothesis values used with
%   the sbox computation are the same as the target values used when
%   creating the mean vectors. That is, the values in dparams.values
%   represent both the target values for profiling and all the possible key
%   hypothesis.
%
%   This method returns two parameters:
%   - d: a vector of discriminant scores (may be probabilities) d(k | X)
%   - p: a vector of probabilities p(X | k), which makes sense mostly if
%   you are using dtype = mvnpdf or dmvnpdfnorm or dmvnpdfnormspool
%
%   Author: Marios O. Choudary (osc22@cl.cam.ac.uk)
%
%   See the "Applied Multivariate Statistical Analysis" book for details.
%
%   See also evaluate_discriminant, evaluate_mdlike, mvnlikelog, mdlike.


%% Check and initialise parameters
% You should use this in the script calling this method: addpath('AES/');
% replacing the path with the place where your AES implementation is.
[s_box, ~] = s_box_gen;
if ~isfield(dparams, 'mvec')
    error('Missing mvec parameter');
end
N = size(X, 1);
D = size(X, 2);
nr_groups = size(dparams.mvec, 1);
if size(dparams.mvec, 2) ~= D
    error('Incompatible size of mvec with X');
end
if length(V) ~= N
    error('Incorrect length of vector V');
end
d = zeros(nr_groups, 1);
p = zeros(nr_groups, 1);

%% Compute scores from given parameters
if strcmp(dtype, 'dlinear')
    %% Compute linear discriminant using common covariance
    ct = -1/2;
    for k = 1:nr_groups
        d(k) = 0;
        for j = 1:N
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            xm = dparams.mvec(idx, :);
            d(k) = d(k) + xm*dparams.sinv*X(j,:)' + ...
                ct*(xm*dparams.sinv*xm');
        
        end        
    end
elseif strcmp(dtype, 'dlinearnocov')
    %% Compute linear discriminant using identity covariance
    ct = -1/2;
    for k = 1:nr_groups
        d(k) = 0;
        for j = 1:N
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            xm = dparams.mvec(idx, :);
            d(k) = d(k) + xm*X(j,:)' + ...
                ct*(xm*xm');        
        end      
    end
elseif strcmp(dtype, 'dlog')
    %% Compute log-likelihood discriminant
    ct = -1/2;
    for k = 1:nr_groups
        dsum = 0;
        for j = 1:N
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            x = X(j,:) - dparams.mvec(idx,:);
            dsum = dsum + dparams.slogdet(idx) + x*dparams.sinv(:,:,idx)*x';
        end
        d(k) = ct * dsum;
    end
elseif strcmp(dtype, 'dmvnpdf')
    %% Compute multivariate normal pdf, going first in log domain
    ct1 = -(D*N/2)*log(2*pi);
    ct2 = -1/2;
    for k = 1:nr_groups
        dsum = 0;
        for j = 1:N
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            x = X(j,:) - dparams.mvec(idx,:);
            dsum = dsum + dparams.slogdet(idx) + x*dparams.sinv(:,:,idx)*x';
        end
        vlog = ct1 + ct2*dsum;
        d(k) = exp(vlog);
    end  
elseif strcmp(dtype, 'dmvnpdfnorm') || strcmp(dtype, 'dmvnpdfnormspool')
    %% Compute mvn pdf, normalising all probabilities after each attack trace
    for k=1:nr_groups
        d(k) = 1;
        p(k) = 1;
    end
    
    for j = 1:N
        for k = 1:nr_groups
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            pmvn = mvnpdfinv(X(j,:), dparams.mvec(idx,:), ...
                dparams.sinv(:,:,idx), dparams.slogdet(idx));
            d(k) = d(k) * pmvn;
            p(k) = p(k) * pmvn;
        end
        sdk = sum(d);        
        for k = 1:nr_groups
            d(k) = d(k) / sdk;
        end
    end    
elseif strcmp(dtype, 'dmdpdf')
    if size(dparams.sinv, 3) > 1
        error('dmdpdf only supported for single pooled covariance');
    end
    %% Compute Mahalanobis distance pdf based on Chi-square PDF
    % Compute first in log domain to reduce numerical errors
    for k = 1:nr_groups
        vlog = 0;
        for j = 1:N
            v = s_box(bitxor(dparams.values(k), V(j))+1);
            idx = (dparams.values == v);
            x = X(j,:) - dparams.mvec(idx,:);
            md = x*dparams.sinv*x';
            vlog = vlog + log(chi2pdf(md, D));
        end
        d(k) = exp(vlog);
    end    
end

end



