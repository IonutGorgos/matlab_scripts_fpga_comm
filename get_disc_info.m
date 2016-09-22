function [disc_info] = get_disc_info(...
                                        X_attack, V_attack, ...
                                        nr_traces_vec, ...
                                        dtype, dparams)
%GET_DISC_INFO Returns discriminant values from generic data
%   [disc_info] = GET_DISC_INFO(...
%                                   X_attack, V_attack, ...
%                                   nr_traces_vec, ...
%                                   dtype, dparams)
%   computes values of some discriminant for generic data.
%
%   X_attack should be a matrix of size nr_traces x nr_samples x nr_groups
%   having the test data which has been already
%   preprocessed (e.g. compressed) to the same format (size) as the
%   parameters sent via dparams (see below).
%
%   V_attack should be a vector of length nr_groups, providing
%   the values (e.g. key byte) corresponding to all traces per group in
%   X_attack.
%
%   nr_traces_vec: the groups of traces to test. nr_traces_vec
%   should be a vector of integers, where each element represents the
%   size of one test group. 
%   E.g. if nr_traces_vec = [5, 10] then 2 tests will be
%   performed. The first test will use 5 random samples and the second
%   test will use 10 random samples.
%
%   dtype specifies the type of discriminant to be used. See compute_scores
%   for the available discriminants.
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
%   - 'sbox_plaintext': used to provide the corresponding plaintext of
%   each attack trace. If given then this should be a matrix of size
%   nr_traces x nr_groups. In this case, this method will use each
%   plaintext to compute the value v=aes-sbox(plaintext XOR key_hypothesis)
%   for each trace and then use the value v in the computation of the
%   discriminants.
%
%   disc_info is a structure containing the results of this
%   method. disc_info has these substructures:
%   - scores: matrix of size nr_values x nr_attack_groups,
%   containing the discriminant scores for all values of the given
%   parameters (e.g. means) and all attack groups. This is d(k | X).
%   - probs: matrix of size nr_values x nr_attack_groups,
%   containing the probabilities for the leakage traces.
%   This is p(X | k). Only available for some cases for now.
%   - depth: vectors of length nr_attack_groups, containing the index
%   of the correct value/group when sorting the scores.
%   This data can be used to compute the guessing entropy for example.
%   - rindex: vectors of length nr_traces, containing the indeces
%   of the random selection for each group size.
%   
%   All structures (e.g. disc_info.scores) have
%   substructures corresponding to the nr_traces_vec
%   parameters. These structures are labeled group1, group2, etc...
%   For example, if nr_traces_vec = [1 10 50] then disc_info.scores
%   will have the following substructures:
%       - group1: results for tests with nr_traces_vec(1) = 1 trace
%       - group2: results for tests with nr_traces_vec(2) = 10 traces
%       - group3: results for tests with nr_traces_vec(3) = 50 traces
%
%   These last structures in turn contain the actual data.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

%% Initialize and check parameters
if ~isfield(dparams, 'mvec')
    error('Missing mvec parameter');
end
nr_traces = size(X_attack, 1);
nr_attack_groups = size(X_attack, 3);
nr_test_groups = length(nr_traces_vec);
nr_values = size(dparams.mvec, 1);

%% Compute the scores for each test group size
for i=1:nr_test_groups
    fprintf('Computing scores for group size %d\n', ...
             nr_traces_vec(i));
    
    %% Set up the score and proba matrices
    disc_info.scores.(['group' num2str(i)]) = ...
        zeros(nr_values, nr_attack_groups);
    disc_info.probs.(['group' num2str(i)]) = ...
        zeros(nr_values, nr_attack_groups);
    
    %% Select traces to use for attack
    if isfield(dparams, 'random')
        rindex = randi([1, nr_traces], nr_traces_vec(i), 1);
    else
        rindex = (1:nr_traces_vec(i))';
    end
    disc_info.rindex.(['group' num2str(i)]) = rindex;
    
    %% Compute scores for each attack group
    for group=1:nr_attack_groups        
        %% Select data
        data = X_attack(rindex, :, group);

        %% Compute scores/probabilities
        if isfield(dparams, 'sbox_plaintext') && ~isempty(dparams.sbox_plaintext)
            plaintext = dparams.sbox_plaintext(rindex,group);
            [d, p] = compute_scores_sbox(data, plaintext, dtype, dparams);
            disc_info.probs.(['group' num2str(i)])(:, group) = p(:);
        else
            d = compute_scores(data, dtype, dparams);
        end
        disc_info.scores.(['group' num2str(i)])(:, group) = d(:);

        %% Then compute depth vectors
        % The main assumption here is that the group corresponding to
        % the highest score is the correct group.
        [~, si] = sort(d(:), 1, 'descend');
        disc_info.depth.(['group' num2str(i)])(group) = ...
            find(dparams.values(si) == V_attack(group));
    end
    toc
end

end
    
    