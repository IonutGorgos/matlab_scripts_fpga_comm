function d = compute_scores_cpa_sbox(X, P, atype, dparams)
%COMPUTE_SCORES_CPA_SBOX Computes CPA scores in an AES-SBOX scenario
%   [d] = COMPUTE_SCORES_CPA_SBOX(X, P, atype, dparams)
%   computes some CPA scores of traces X based on the given parameters.
%
%   X should be a matrix of size N x D, having N D-dimensional observations
%   for which the scores will be computed.
%
%   P should be a vector of length N (number of traces), where each element
%   is the plaintext corresponding to a trace. For each key hypothesis k*,
%   a value v=aes-sbox(plaintext XOR k*) will be computed and this value
%   will be used for each trace in the computation of the discriminants.
%
%   atype specifies the type of CPA attack/score to be used. Currently defined
%   are:
%   - 'classic': the standard CPA approach (see Brier et al. '04) with normalisation
%   - 'zcdf': using a z-transform and then a pdf/cdf
%
%   dparams is a structure that provides the necessary parameters for
%   computing the scores. These paramters may include:
%   - 'values': a vector of length nr_values specifying the key hypothesis
%   values to be tested during the attack.
%   - 'lmodel': leakage model for each possible S-box output value. This
%   should basically be a vector of length 256.
%
%   This method returns a matrix d of size 256 x D, having the scores (some
%   sort of normalised correlation) for each possible S-box output value
%   and leakage sample.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)
%
%   See also evaluate_discriminant, evaluate_mdlike, mvnlikelog, mdlike.


%% Check and initialise parameters
addpath('AES/');
[s_box, ~] = s_box_gen;
if ~isfield(dparams, 'values')
    error('Missing values parameter');
end
nr_values = length(dparams.values);
lmodel = dparams.lmodel;
N = size(X, 1);
D = size(X, 2);
if length(P) ~= N
    error('Incorrect length of vector V');
end
d = zeros(nr_values, D);

%% Compute classic CPA from given parameters
if strcmp(atype, 'classic')
    %% Compute classic CPA with normalisation of scores
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            d(k, i) = abs(c(1,2));
        end
        % Normalise correlation
        dsum = sum(d(:,i));
        d(:,i) = d(:,i) ./ dsum;
    end
elseif strcmp(atype, 'zcdf')
    %% Compute CPA scores via z-transform and then CDF
    zvar = 1 / (N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = abs(c(1,2));
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
    end
elseif strcmp(atype, 'zcdfnorm')
    %% Compute CPA scores via z-transform and then CDF, then normalize
    zvar = 1 / (N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = abs(c(1,2));
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
        % Normalise correlation (a CDF does not imply normalisation)
        dsum = sum(d(:,i));
        d(:,i) = d(:,i) ./ dsum;
    end
elseif strcmp(atype, 'zcdfgood')
    %% Compute CPA scores via z-transform with good params and then CDF
    zvar = 1 / (N-3);
    zvec = zeros(nr_values, 1);
    for i=1:D
        %% First compute the z values
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = abs(c(1,2));
            zvec(k) = 0.5 * log((1+c)/(1-c));            
        end
        %% Compute mean of Z values
        zmean = mean(zvec);
        %% Now compute the transformed values
        for k = 1:nr_values
            d(k, i) = mvncdf(zvec(k), zmean, zvar);
        end
        %% Normalise correlation (a CDF does not imply normalisation)
        dsum = sum(d(:,i));
        d(:,i) = d(:,i) ./ dsum;
    end
elseif strcmp(atype, 'zcdfbest')
    %% Compute CPA scores via z-transform with best params and then CDF
    zvar = 1;
    zvec = zeros(nr_values, 1);
    for i=1:D
        %% First compute the z values
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = abs(c(1,2));
            zvec(k) = 0.5 * log((1+c)/(1-c));            
        end
        %% Compute mean of Z values
        zmean = mean(zvec);
        %% Now compute the transformed values
        for k = 1:nr_values
            d(k, i) = mvncdf(zvec(k), zmean, zvar);
        end
        %% Normalise correlation (a CDF does not imply normalisation)
        dsum = sum(d(:,i));
        d(:,i) = d(:,i) ./ dsum;
    end
elseif strcmp(atype, 'zcdfsqrt')
    %% Compute CPA scores via z-transform and then CDF
    zvar = 1 / sqrt(N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = abs(c(1,2));
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
    end
elseif strcmp(atype, 'zcdfnoabs')
    %% Compute CPA scores via z-transform of abs value and then CDF
    zvar = 1 / (N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = c(1,2);
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
    end
elseif strcmp(atype, 'zcdfnoabsnorm')
    %% Compute CPA scores via z-transform of abs value and then CDF
    zvar = 1 / (N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = c(1,2);
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
        % Normalise correlation (a CDF does not imply normalisation)
        dsum = sum(d(:,i));
        d(:,i) = d(:,i) ./ dsum;
    end
elseif strcmp(atype, 'zcdfnoabssqrt')
    %% Compute CPA scores via z-transform of abs value and then CDF
    zvar = 1 / sqrt(N-3);
    for i=1:D
        for k = 1:nr_values
            % Need (+1) in the 2 lines below due to Matlab indexing from 1
            V = s_box(bitxor(dparams.values(k), P) + 1);
            L = lmodel(V + 1);
            c = corrcoef(X(:,i), L);
            c = c(1,2);
            z = 0.5 * log((1+c)/(1-c));
            d(k, i) = mvncdf(z, 0, zvar);
        end
    end    
end

end



