% Test profiled attacks on AES
% Author: Omar Choudary

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none'); % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
addpath('/home/osc22/projects/power_analysis_xmega/mscripts/');
addpath('/home/osc22/projects/power_analysis_xmega/mscripts/AES/');
data_path = 'polar/';
fmap_profile = [data_path, 'results/raw/p2/beta/p2_bat_fb_2mhz_beta_k0.raw'];
fmap_attack = [data_path, 'results/raw/p2/beta/p2_bat_fb_2mhz_beta_k1.raw'];
data_title = 'Templates P2 BAT FB';
path_data = [data_path, 'results/a2/'];
name_data = sprintf('templates_aes_p2_beta_b1_dlinear_all_n1000_g1000_r10.mat');
n_profile = 1000;
nr_traces_vec = [1:10, 20:10:100, 200:100:1000];
nr_iterations = 10;
samples = 1:5000;
rng('default'); % use same randomisation to get consistent results

%% Load files
fprintf('Mapping data\n');
[m_data_profile, metadata_profile] = get_mmap(fmap_profile);
[m_data_attack, metadata_attack] = get_mmap(fmap_attack);
toc

%% Select target byte
target_byte = 1;

%% Extract profile data
fprintf('Obtaining profile data...\n');
P = double(m_data_profile.data(1).P(target_byte,:)');
key_profile = double(m_data_profile.data(1).K(target_byte));
X_profile = double(m_data_profile.data(1).X(samples,:)');

%% Select combination method
comb_func = @combine_aes_sbox;

%% Pre-compute the S-box for faster computations
[s_box, ~] = s_box_gen;
cf_params = [];
cf_params.s_box = s_box;

%% Compute target values for profiling traces
X_profile_target = comb_func(P, key_profile, cf_params);

%% Extract values for profile/discriminant
V_profile = unique(X_profile_target);
nr_values = length(V_profile);
nr_samples = size(X_profile, 2);

%% Select target bytes (keys)
% Targetting only first AES key byte
key_attack = double(m_data_attack.data(1).K(target_byte));
V_attack = key_attack;

%% Select attack data
fprintf('Obtaining attack data...\n');
X_attack_input = double(m_data_attack.data(1).P(target_byte,:)');
X_attack = double(m_data_attack.data(1).X(samples,:)');

%% Select possible key values
V_key = 0:255;

%% Run attack for each iteration and experiment
nr_attack_groups = length(nr_traces_vec);
nr_exp = 6; % LDA, PCA, 1ppc, 3ppc, 20ppc, allap
results = cell(nr_iterations, nr_exp);
for i=1:nr_iterations
    fprintf('Running attack for iteration i=%d\n', i);
    
    %% Select profile traces
    x_profile_grouped = zeros(n_profile, nr_samples, nr_values);
    for k=1:nr_values
        idx = find(X_profile_target == V_profile(k));
        total_idx = length(idx);
        idx_profile = randperm(total_idx, n_profile);        
        idx = idx(idx_profile);
        x_profile_grouped(:,:,k) = X_profile(idx, :);
    end   
    
    %% Run attack
    fprintf('Running attack for LDA...\n');
    cmethod = 'LDA';
    cparams = [];
    cparams.lda_threshold = 0.95;
    cparams.lda_dimensions = 10;
    dtype = 'dlinearnocov';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 1} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc   
    
    %% Run attack for PCA
    fprintf('Running attack for PCA...\n');
    cmethod = 'PCA';
    cparams = [];
    cparams.pca_threshold = 0.95;
    cparams.pca_alternate = 0;
    cparams.pca_dimensions = 10;
    dtype = 'dlinear';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 2} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc  
    
    %% Run attack for selection
    fprintf('Running attack for 1ppc...\n');
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'ftest';
    cparams.sel = '1ppc';
    cparams.p1 = 230;
    dtype = 'dlinear';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 3} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc  
    
     %% Run attack for selection
    fprintf('Running attack for 3ppc...\n');
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'ftest';
    cparams.sel = '3ppc';
    cparams.p1 = 230;
    dtype = 'dlinear';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 4} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc  
    
    %% Run attack for selection
    fprintf('Running attack for 20ppc...\n');
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'ftest';
    cparams.sel = '20ppc';
    cparams.p1 = 230;
    dtype = 'dlinear';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 5} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc  
    
    %% Run attack for selection
    fprintf('Running attack for allap...\n');
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'ftest';
    cparams.sel = 'allap';
    cparams.p1 = 0.95;
    dtype = 'dlinear';
    eparams = [];
    eparams.cf_params = cf_params;
    eparams.save_eval = 1;
    results{i, 6} = run_template_attack_comb(...
        x_profile_grouped, V_profile, ...
        X_attack, X_attack_input, ...
        V_attack, V_key, comb_func, ...
        cmethod, cparams, ...
        dtype, ...
        nr_traces_vec, eparams);
    toc  
end

%% Save all variables and clean up
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results', '-v7.3');
toc

%% Exit when running in script mode
% exit
