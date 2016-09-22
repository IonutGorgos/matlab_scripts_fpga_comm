%% Test templates
% Author: Omar Choudary

%% Reset environment
close all;
clear;
set(0, 'DefaulttextInterpreter', 'none') % Remove TeX interpretation
tic

%% Setup the necessary paths and parameters
addpath('/home/osc22/projects/power_analysis_xmega/mscripts/');
% data_path = '/anfs/bigdisc/osc22/data/grizzly/';
data_path = '/media/duba/data/grizzly/';
fmap_alpha = [data_path, 'results/raw/e2/alpha/e2_bat_fb_alpha_raw_s_0_3071.raw'];
fmap_beta = [data_path, 'results/raw/e2/beta/e2_bat_fb_beta_raw_s_0_3071.raw'];
fmap_gamma = [data_path, 'results/raw/e2/gamma/e2_bat_fb_gamma_raw_s_0_3071.raw'];
fmap_delta = [data_path, 'results/raw/e2/delta/e2_bat_fb_delta_raw_s_0_3071.raw'];
data_title = 'Templates A2D';
path_data = [data_path, 'results/test/'];
name_data = sprintf('a2d_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat');
rand_iter = 1; %0;
n_profile = 200; % ensure that: n_profile + n_attack < nr_blocks
nr_traces_vec = [1:10]; %, 20:10:100, 200, 500, 1000];
bytes = 0:255;
atype = 'mvn';

%% Load files for all boards
fprintf('Mapping data for Alpha\n');
[m_data_alpha, metadata_alpha] = get_mmap(fmap_alpha);
toc

fprintf('Mapping data for Beta\n');
[m_data_beta, metadata_beta] = get_mmap(fmap_beta);
toc

fprintf('Mapping data for Gamma\n');
[m_data_gamma, metadata_gamma] = get_mmap(fmap_gamma);
toc

fprintf('Mapping data for Delta\n');
[m_data_delta, metadata_delta] = get_mmap(fmap_delta);
toc

%% Put data into cells
mmap_data = cell(4,1);
metadata = cell(4,1);
names = cell(4,1);
mmap_data{1} = m_data_alpha;
metadata{1} = metadata_alpha;
names{1} = 'Alpha';
mmap_data{2} = m_data_beta;
metadata{2} = metadata_beta;
names{2} = 'Beta';
mmap_data{3} = m_data_gamma;
metadata{3} = metadata_gamma;
names{3} = 'Gamma';
mmap_data{4} = m_data_delta;
metadata{4} = metadata_delta;
names{4} = 'Delta';

%% Select idx for profile/attack
nr_blocks = 3072;
idx = 1:nr_blocks;
idx_profile = union(find(mod(idx,3) == 1), find(mod(idx,3) == 2));
idx_profile = idx_profile(randi([1, length(idx_profile)], 1, n_profile));
idx_attack = find(mod(idx,3) == 0);

%% Setup attack combinations of devices
% profile_dev = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};
% attack_dev = {2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3};
profile_dev = {1};
attack_dev = {2};
nr_exp = length(profile_dev);

%% Set up attack/result cells
results = cell(nr_exp, 6);

%% Run attack for each combination of profile/attack data
for k=1:nr_exp
    
    %% Run attack for LDA, K=4
    cmethod = 'LDA';
    cparams = [];
    cparams.lda_dimensions = 4;
    cparams.lda_threshold = 0.95;
    discriminant = 'linearnocov';
    eparams = [];    
    results{k,1} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);

    %% Run attack for PCA, K=4
    cmethod = 'PCA';
    cparams = [];
    cparams.pca_threshold = 0.95;
    cparams.pca_alternate = 0;
    cparams.pca_dimensions = 4;
    discriminant = 'linear';
    eparams = [];
    results{k,2} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);

    %% Run attack for 1ppc, using adaptation
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'dom';
    cparams.sel = '1ppc';
    cparams.p1 = 240;
    discriminant = 'linear';
    eparams = [];
    results{k,3} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);

    %% Run attack for 3ppc
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'dom';
    cparams.sel = '3ppc';
    cparams.p1 = 240;
    discriminant = 'linear';
    eparams = [];
    results{k,4} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);

    %% Run attack for 20ppc
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'dom';
    cparams.sel = '20ppc';
    cparams.p1 = 240;
    discriminant = 'linear';
    eparams = [];
    results{k,5} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);

    %% Run attack for allap
    cmethod = 'sample';
    cparams = [];
    cparams.curve = 'dom';
    cparams.sel = 'allap';
    cparams.p1 = 0.95;
    discriminant = 'linear';
    eparams = [];
    results{k,6} = run_template_attack(...
        mmap_data{profile_dev{k}}, metadata{profile_dev{k}}, idx_profile, ...
        mmap_data{attack_dev{k}}, metadata{attack_dev{k}}, idx_attack, ...
        bytes, atype, cmethod, cparams, discriminant, ...
        rand_iter, nr_traces_vec, eparams);
end
    
%% Save results
fprintf('All done, saving data...\n');
save([path_data, name_data], 'results', 'names', 'profile_dev', 'attack_dev', '-v7.3');
toc

%% Exit, only use for condor. Matlab will quit.
% exit

