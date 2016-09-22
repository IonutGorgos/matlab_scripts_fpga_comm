% Test key enumeration
% Author: Omar Choudary

%% Reset environment
%close all;
%clear;
set(0, 'DefaulttextInterpreter', 'none'); % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
% addpath('mscripts/');
addpath('key_enum/'); % path for keyEnumeration and generateLeakageMatlab
addpath('AES/');
%data_title = 'Scores SIM DATA';
% data_path = 'grizzly/';
% path_results = [data_path, 'results/keyenum/'];
%data_path = 'key_enum/';
path_results = [data_path, 'results/'];
%name_results = sprintf('keyenum_gendata_cpa_classic_var1_g50_r50.mat');
path_data = ['key_enum/data/'];
%name_data = sprintf('simdata_var1_np50000_na5000_m1_r100_v2.mat');
nr_attack_bytes = 4;
nr_traces_vec = [20:5:100]; % Use this for var=10
% nr_traces_vec = [10:20, 30:10:50]; % Use this for var=1
nr_iterations = 50;
rng('default'); % use same randomisation to get consistent results

%% Set possible candidate values
target_values = 0:255;
nr_values = length(target_values);

%% Set Hamming Weight as leakage model for each value in simulated data
lmodel = hamming_weight(target_values);

%% Load previously generated data
% For each possible attack iteration, simdata{i} contains the following:
% 'M_profile': vector of plaintexts for profiling
% 'X_profile' vector of leakage traces for profiling
% 'K_profile': key used during profiling (just 1 byte)
% 'M_attack': vector of plaintexts for attack
% 'X_attack': matrix of leakage for attack, for all 4 bytes (see below)
% 'K_attack': the four key byte values for the attack data
%load([path_data, name_data]);

%% Generate sbox
[s_box, ~] = s_box_gen;

%% Run key enumeration for each iteration
nr_attack_groups = length(nr_traces_vec);
keypos = zeros(nr_attack_groups, nr_iterations);
results = cell(nr_iterations, nr_attack_bytes);
for i=1:nr_iterations
    fprintf('Running key enumeration on sim data for i=%d\n', i);
    
    %% Run attack for each target byte
    for k=1:nr_attack_bytes
        fprintf('Running attack for target byte %d...\n', k);
        atype = 'classic';
        eparams = [];
        results{i,k} = run_cpa_sbox_scores(...
            simdata{i}.X_attack(:,:,k), simdata{i}.M_attack, ...
            simdata{i}.K_attack(k), ...
            lmodel, atype, nr_traces_vec, eparams);
        toc
    end
    
end

%% Save all variables and clean up
fprintf('All done, saving data...\n');
save([path_results, name_results], 'results', 'keypos', '-v7.3');
toc

%% Exit when running in script mode
% exit
