function [results] = run_template_attack_sbox_scores(...
    X_profile, X_profile_target, ...
    X_attack, X_attack_input, ...
    V_profile, V_attack, ...
    dtype, ...
    nr_traces_vec, eparams)
%RUN_TEMPLATE_ATTACK_SBOX_SCORES Runs a template attack and returns scores
%   [results] = RUN_TEMPLATE_ATTACK_SBOX_SCORES(...
%       X_profile, X_profile_target, ...
%       X_attack, X_attack_input, ...
%       V_profile, V_attack, ...
%       dtype, ...
%       nr_traces_vec, eparams)
%   runs a template attack on data given as matrices for profile and attack
%   and returns a results scores/probabilities of each possible value/group.
%
%   X_profile should be a matrix of size nr_trials x nr_samples, where
%   nr_trials is the number of traces each having nr_samples leakage
%   samples. This matrix should contain the data for profiling.
%
%   X_profile_target should be a vector of length nr_trials, containing the
%   target value (this method assumes target=S-box(key XOR plaintext)) for
%   each leakage trace, so that the profiling can be done for each target
%   value by grouping traces.
%
%   X_attack should be a matrix of size nr_trials x nr_samples x nr_bytes
%   containing the leakage traces for attacking different key byte values
%   (one nr_trials x nr_samples matrix per key byte value). This method
%   will compute the discriminants/scores per key byte value.
%
%   X_attack_input should be a matrix of size nr_trials x nr_bytes
%   specifying the plaintext used for each of the leakage traces in the
%   X_attack matrix. During the template attack, this method will compute
%   the value S-box(key XOR plaintext) for each key byte hypopthesis.
%
%   V_profile should be a vector of length nr_values, specifying the
%   (target) values for which profiling should be performed. These values
%   should be found within X_profile_target.
%
%   V_attack should be a vector of length nr_bytes specifying the
%   actual values of the key bytes for which the matrix X_attack was given.
%
%   dtype specifies the type of discriminant to be used. See compute_scores
%   for available types.
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used for each element.
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
nr_samples = size(X_profile, 2);
nr_values = length(V_profile);
results = [];
results.dtype = dtype;
results.nr_traces_vec = nr_traces_vec;
results.v_attack = V_attack;
if nargin < 9
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_template_attack_sbox_scores() ...\n');


%% Compute templates and evaluation parameters
fprintf('Computing template and evaluation parameters...\n');
tmiu = zeros(nr_values, nr_samples);
tsigma = zeros(nr_samples, nr_samples, nr_values);
for k=1:nr_values
    idx = X_profile_target == V_profile(k);
    x = X_profile(idx, :);
    tmiu(k,:) = mean(x, 1);
    tsigma(:,:,k) = cov(x);
end

dparams = [];
dparams.values = V_profile;
dparams.mvec = tmiu;
dparams.scov = tsigma;
dparams.sbox_plaintext = X_attack_input;

% TODO: implement/add for linearfast
if (strcmp(dtype, 'dlinear') || strcmp(dtype, 'dmdpdf'))
    c0 = mean(tsigma, 3);
    ic0 = inv(c0);
    dparams.sinv = ic0;
elseif (strcmp(dtype, 'dlog') || strcmp(dtype, 'dmvnpdf') || strcmp(dtype,'dmvnpdfnorm'))
    n = size(tsigma, 3);
    tsinv = zeros(size(tsigma));
    tlogdet = zeros(n, 1);
    for k = 1:n
        tsinv(:,:,k) = inv(tsigma(:,:,k));
        tlogdet(k) = logdet(tsigma(:,:,k), 'chol');
    end
    dparams.sinv = tsinv;
    dparams.slogdet = tlogdet;
elseif strcmp(dtype, 'dmvnpdfnormspool')
    c0 = mean(tsigma, 3);
    ic0 = inv(c0);
    l0 = logdet(c0, 'chol');
    n = size(tsigma, 3);
    tsinv = repmat(ic0, 1, 1, n);
    tlogdet = repmat(l0, n, 1);
    dparams.sinv = tsinv;
    dparams.slogdet = tlogdet;    
end

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.dparams = dparams;
end

%% Compute the discriminant scores/probabilities
fprintf('Computing discriminant info...\n');
[results.disc_info] = get_disc_info(...
                               X_attack, V_attack,...
                               nr_traces_vec, ...
                               dtype, dparams);
toc

end
