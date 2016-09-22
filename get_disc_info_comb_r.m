function [disc_info] = get_disc_info_comb_r(...
                                        X_attack, V_attack, V_key, ...
                                        nr_traces_vec, nr_iter, ...
                                        dtype, dparams, comb_func)
%GET_DISC_INFO Returns discriminant values from generic data
%   [disc_info] = GET_DISC_INFO_COMB(...
%                                   X_attack, V_attack, ...
%                                   nr_traces_vec, ...
%                                   dtype, dparams, comb_func)
%   computes values of a discriminant for some target data.
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
%   V_key should be a vector of length nr_keys specifying all the possible
%   value of the key, i.e. this should be the vector of key hypothesis,
%   that will be used to find the most likely value.
%
%   nr_traces_vec: the groups of traces to test. nr_traces_vec
%   should be a vector of integers, where each element represents the
%   size of one test group. 
%   E.g. if nr_traces_vec = [5, 10] then 2 tests will be
%   performed. The first test will use 5 random samples and the second
%   test will use 10 random samples.
%
%   nr_iter specifies how many times to run each test on
%   randomly picked samples to produce the desired results.
%   Note that the execution time increases linearly with this parameter.
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
%   - 'plaintext': used to provide the corresponding plaintext of
%   each attack trace. If given then this should be a matrix of size
%   nr_traces x nr_groups. This method will use each
%   plaintext to compute the value v=comb_func(plaintext, key_hypothesis)
%   for each trace and then use the value v in the computation of the
%   discriminants.
%   - 'cf_params': parameters for the combination function.
%
%   comb_func should be a method, taking two parameters: a vector of
%   plaintext values and a vector of key bytes, and which should return the
%   target combination of these (e.g. Sbox(p XOR k)). This combination
%   function should expect values starting with 0.
%
%   disc_info is a structure containing the results of this
%   method. disc_info has these substructures:
%   - scores: matrix of size nr_keys x nr_attack_groups,
%   containing the discriminant scores for all keys and all attack groups.
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
nr_keys = length(V_key);

%% Compute the scores for each test group size
for i=1:nr_test_groups
    fprintf('Computing scores for group size %d\n', ...
             nr_traces_vec(i));
    
    %% Set up the score matrices
    disc_info.scores.(['group' num2str(i)]) = ...
        zeros(nr_keys, nr_attack_groups, nr_iter);
    
    %% Set up the depth vectors
    disc_info.depth.(['group' num2str(i)]) = zeros(nr_attack_groups, nr_iter);
    
    %% Select random traces to use for all the tests
    rindex = zeros(nr_traces_vec(i), nr_iter);
    for count=1:nr_iter    
        rindex(:,count) = randperm(nr_traces, nr_traces_vec(i))';
    end
    disc_info.rindex.(['group' num2str(i)]) = rindex;

    %% Compute scores for each attack group
    for group=1:nr_attack_groups   
        %% Perform the tests for nr_iter
         for count=1:nr_iter   
           %% Select data
            data = X_attack(rindex(:,count), :, group);

           %% Compute scores/probabilities        
            plaintext = dparams.plaintext(rindex(:,count),group);
            d = compute_scores_comb(data, plaintext, V_key, dtype, dparams, comb_func);
            disc_info.scores.(['group' num2str(i)])(:, group, count) = d(:);

           %% Then compute depth vectors
            % The main assumption here is that the group corresponding to
            % the highest score is the correct group.
            [~, si] = sort(d(:), 1, 'descend');
            disc_info.depth.(['group' num2str(i)])(group,count) = ...
                find(V_key(si) == V_attack(group));
        end
    end
    toc
end

end
    
    