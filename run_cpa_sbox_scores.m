function [results] = run_cpa_sbox_scores(...
    X_attack, X_attack_input, key_val, ...
    lmodel, atype, nr_traces_vec, eparams)
%RUN_CPA_SBOX_SCORES Runs a CPA attack and returns scores for candidates
%   [results] = RUN_CPA_SBOX_SCORES(...
%       X_attack, X_attack_input, key_val, ...
%       lmodel, atype, nr_traces_vec, eparams)
%   runs a CPA attack on data given as matrices
%   and returns a results scores/probabilities of each possible value/group.
%
%   X_attack should be a matrix of size nr_trials x nr_samples
%   containing the leakage traces for attack.
%
%   X_attack_input should be a vector of length nr_trials
%   specifying the plaintext used for each of the leakage traces in the
%   X_attack matrix. During the attack, this method will compute
%   the value S-box(key XOR plaintext) for each key byte hypopthesis.
%
%   key_val should specify the actual value of the key byte for which
%   the matrix X_attack was given.
%
%   lmodel should be a vector of length 256, representing the leakage for
%   each possible S-box output value.
%
%   atype specifies the type of attack performed. See the script
%   'compute_scores_cpa_sbox' for the complete list of possible types.
%   Some of the supported types are:
%   - 'classic': the standard CPA approach (see Brier et al. '04)
%   - 'zcdf': using a z-transform and then a cdf
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used during the attacks.
%
%   eparams is a structure of extra parameters that may be needed for some
%   options.
%
%   This function returns a structure 'results' with the following:
%   -> results.atype: the attack type string.
%   -> results.nr_traces_vec: the vector with number of attack traces per
%   attack. That is, we can have several attacks, each with a different
%   number of traces.
%   -> results.disc_info: scores/probabilities for each leakage sample and
%      candidate key value. This is a matrix of size
%      256 x nr_samples x len(nr_traces_vec) for each group of attack
%      traces.
%   -> results.error: an optional messsage in case of an error. This may
%   avoid the abortion of a larger set of runs while allowing to detect the
%   cases in which an error occurred.
%
%   See the paper "Correlation Power Analysis with a Leakage Model", Brier
%   et al., 2004.

%% Check and initialise parameters
addpath('AES/');
nr_traces = size(X_attack, 1);
nr_samples = size(X_attack, 2);
nr_test_groups = length(nr_traces_vec);
results = [];
results.atype = atype;
results.nr_traces_vec = nr_traces_vec;
nr_values = 256; % We know we are dealing with an 8-bit S-box
key_val_candidates = 0:255;
if nargin < 6
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_cpa_sbox_scores() ...\n');

%% Compute the scores for each test group size
for i=1:nr_test_groups
    fprintf('Computing scores for group size %d\n', ...
             nr_traces_vec(i));
    
    %% Set up the score/depth matrices
    results.disc_info.scores.(['group' num2str(i)]) = zeros(nr_values, nr_samples);
    results.disc_info.depth.(['group' num2str(i)]) = zeros(nr_samples);
    
    %% Select traces to use for all the tests
    if isfield(eparams, 'random')
        rindex = randi([1, nr_traces], nr_traces_vec(i), 1);
    else
        rindex = 1:nr_traces_vec(i);
    end
    results.disc_info.rindex.(['group' num2str(i)]) = rindex;
    
    %% Select data
    data = X_attack(rindex,:);
    
    %% Set/specify parameters
    dparams = [];
    dparams.values = key_val_candidates;
    dparams.lmodel = lmodel;

    %% Compute scores/probabilities
    plaintext = X_attack_input(rindex,:);
    d = compute_scores_cpa_sbox(data, plaintext, atype, dparams);
    results.disc_info.scores.(['group' num2str(i)]) = d;

    %% Then compute depth vectors for each leakage sample    
    for idx_sample=1:nr_samples
        % The main assumption here is that the group corresponding to
        % the highest score is the correct group.
        [~, si] = sort(d(:,idx_sample), 1, 'descend');
        results.disc_info.depth.(['group' num2str(i)])(idx_sample) = ...
            find(dparams.values(si) == key_val);
        
    end
    toc
end

end
