%   this script is used to submit mutant and wt jobs to orchestra
%   for mig1d, we force aR=0 and don't let it change;
%   for gal80d, we force a80=ag80=0, and keep them unchange;
%   the other settings are the same to wt

n_init = 100;     % 5 replicates
% n_propose = 1000000;     % run for 1000,000 iterations
n_propose = 50000;

% folder = '../metaData/ReverseFitting/';
folder = '../metaData/fixOverlap/';

%% use intracellular glucose concentration to replace all Mig1* in the equations
base_param = set_parameter(8);  
parameter_update = readtable('Mar3rd_param_config_set8.csv');
% load('../metaData/noneCompetNoneNoisy.mat');
% load('../metaData/noneCompetAddNoise.mat');
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'wildtype_96well', folder);
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'BC187_Kayla', folder);

%% fit all 96well data, mig1d
load('../metaData/trait_extraction/S288C-double_gradient/mig1d_all_data.mat')
% parameter_update = readtable('Sep7th_param_config_set8.csv');
mig1d_param = base_param;
mig1d_param.aR = 0;
% mig1d_param.aR = 0.001;
% parameter_update{4,'prior_mean'} = 0.001;
% parameter_update{5,'prior_sigma'} = 1;    % no constrain on aR
% parameter_update{5,'proposal_sigma_1'} = 0.7;     % in case the initial value is ridiculous, a bigger step size is needed
rand_init_generator(n_init, trait, n_propose, mig1d_param, parameter_update, 'mig1d_96well', folder);

%% fit all 96well data, gal80d
load('../metaData/trait_extraction/S288C-double_gradient/gal80d_all_data.mat')
parameter_update = readtable('Mar3rd_param_config_set8.csv');
% parameter_update = readtable('Sep7th_param_config_set8.csv');
base_param = set_parameter(8);
gal80d_param = base_param;
gal80d_param.a80 = 0;
gal80d_param.ag80 = 0;
% gal80d_param.a80 = 0.001;
% gal80d_param.ag80 = 0.001;
% parameter_update{[3,8],'prior_mean'} = 0.001;
% parameter_update{[4,11], 'prior_sigma'} = 1;    % no constrain on a80, ag80
% parameter_update{[4,11], 'proposal_sigma_1'} = 0.7;
rand_init_generator(n_init, trait, n_propose, gal80d_param, parameter_update, 'gal80d_96well', folder);


%%
% fit one column
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1c.mat')
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1c', folder);
% fit one row
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r.mat')
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r', folder);
% fit one cross
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r1c.mat')
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r1c', folder);



%% change the formula of Mig1, medium step size
base_param = set_parameter(5);
parameter_update = readtable('Aug8th_param_config_medium_set5.csv');
% % fit one column
% load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1c.mat')
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1c', folder);
% % fit one row
% load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r.mat')
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r', folder);
% fit one cross
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r1c.mat')
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r1c', folder);

% fit GAL80 delete strain, as well as GAL3 and GAL4 level
gal80d_param = base_param;
gal80d_param.a80 = 0.0001;
gal80d_param.ag80 = 0.0001;
load('../metaData/trait_extraction/S288C-double_gradient/gal80d_1c.mat')
rand_init_generator(n_init, trait, n_propose, gal80d_param, parameter_update, 'medium-gal80d_1c', folder);
load('../metaData/trait_extraction/S288C-double_gradient/gal80d_1r.mat')
rand_init_generator(n_init, trait, n_propose, gal80d_param, parameter_update, 'medium-gal80d_1r', folder);

% fit Mig1 delete strain, as well as GAL3 and GAL4 level
mig1d_param = base_param;
mig1d_param.aR = 0.0001;
load('../metaData/trait_extraction/S288C-double_gradient/mig1d_1c.mat')
rand_init_generator(n_init, trait, n_propose, mig1d_param, parameter_update, 'medium-mig1d_1c', folder);
load('../metaData/trait_extraction/S288C-double_gradient/mig1d_1r.mat')
rand_init_generator(n_init, trait, n_propose, mig1d_param, parameter_update, 'medium-mig1d_1r', folder);


%% use alpha*KMglu to replace KMgal, test different step size
% also use beta*kglu to replace kgal
base_param = set_parameter(4);

% fit one column
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1c.mat')

% fit GAL1,3,4; let hill coefficients vary; small step size
parameter_update = readtable('Aug5th_param_config_small_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'small-wildtype_1c', folder);

% fit GAL1,3,4; let hill coefficients vary; medium step size
parameter_update = readtable('Aug5th_param_config_medium_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1c', folder);

% fit GAL1,3,4; let hill coefficients vary; large step size
parameter_update = readtable('Aug5th_param_config_large_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'large-wildtype_1c', folder);


% fit one row
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r.mat')

% fit GAL1,3,4; let hill coefficients vary; small step size
parameter_update = readtable('Aug5th_param_config_small_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'small-wildtype_1r', folder);

% fit GAL1,3,4; let hill coefficients vary; medium step size
parameter_update = readtable('Aug5th_param_config_medium_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r', folder);

% fit GAL1,3,4; let hill coefficients vary; large step size
parameter_update = readtable('Aug5th_param_config_large_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'large-wildtype_1r', folder);


% fit one cross
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r1c.mat')

% fit GAL1,3,4; let hill coefficients vary; small step size
parameter_update = readtable('Aug5th_param_config_small_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'small-wildtype_1r1c', folder);

% fit GAL1,3,4; let hill coefficients vary; medium step size
parameter_update = readtable('Aug5th_param_config_medium_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'medium-wildtype_1r1c', folder);

% fit GAL1,3,4; let hill coefficients vary; large step size
parameter_update = readtable('Aug5th_param_config_large_set4.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'large-wildtype_1r1c', folder);


%% generate config .mat files from MAP parameters
if 0
    folder_MAP = '../results/GAL234-goodfittings_170707/';
    folder_random_init = '../metaData/random_init_mutant_and_wt';
    files = dir(fullfile(folder_MAP, '*.mat'));
    n_file = length(files);
    error_tol = .15;
    parameter_update = readtable('July2nd_parameter_config_wt_set1.csv');
    
    for i_file = 1:n_file
        filename = files(i_file).name;
        load(fullfile(folder_MAP, filename))
        n_propose = 1000000;    % load previous result will overwrite the value
        param_init = param_map;
        
        load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1c.mat')  % trait
        jobtag = 'MAP-wildtype_1c';
        save(fullfile(folder_random_init, [jobtag, num2str(i_file, '_%03d'), '.mat'])...
            , 'parameter_update', 'param_init' ...
            , 'error_tol', 'n_propose' ...
            , 'trait' ...
            , 'jobtag'...
            )
        
        load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r.mat')
        jobtag = 'MAP-wildtype_1r';
        save(fullfile(folder_random_init, [jobtag, num2str(i_file, '_%03d'), '.mat'])...
            , 'parameter_update', 'param_init' ...
            , 'error_tol', 'n_propose' ...
            , 'trait' ...
            , 'jobtag'...
            )
        
        load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r1c.mat')
        jobtag = 'MAP-wildtype_1r1c';
        save(fullfile(folder_random_init, [jobtag, num2str(i_file, '_%03d'), '.mat'])...
            , 'parameter_update', 'param_init' ...
            , 'error_tol', 'n_propose' ...
            , 'trait' ...
            , 'jobtag'...
            )
        
    end
end

%% generate config .mat files that fit BC&YJM single gradient data
base_param = set_parameter(1);
parameter_update = readtable('July2nd_parameter_config_wt_set1.csv');
folder = '../metaData/BCandYJM_single_grad/';

load('../metaData/trait_extraction/BC&YJM-gluc_gradient/BC187_rep_01.mat')
jobtag = 'BC187';
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, jobtag, folder);

load('../metaData/trait_extraction/BC&YJM-gluc_gradient/YJM978_rep_01.mat')
jobtag = 'YJM978';
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, jobtag, folder);


%% generate config .mat files that only fit one column

% load('../metaData/trait_extraction/mig1d_1c.mat')
% base_param = set_parameter(1);
% base_param.aR = 0;
% parameter_update = readtable('MCMC_parameter_config_mig1d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'mig1d_1c');
%
% load('../metaData/trait_extraction/gal80d_1c.mat')
% base_param = set_parameter(1);
% base_param.a80 = 0;
% base_param.ag80 = 0;
% parameter_update = readtable('MCMC_parameter_config_gal80d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'gal80d_1c');

load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1c.mat')
% parameter_update = readtable('July2nd_parameter_config_wt_set1.csv');

% fit GAL1,3,4; let hill coefficients vary
base_param = set_parameter(1);
parameter_update = readtable('Aug3rd_parameter_config_wt_set1.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'varyN-wildtype_1c', folder1);

% set all hill coefficients equal to 1
base_param = set_parameter(3);
parameter_update = readtable('Aug3rd_parameter_config_wt_set3.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'sequestrate-wildtype_1c', folder2);

%% generate config .mat files that only fit one row

% load('../metaData/trait_extraction/mig1d_1r.mat')
% base_param = set_parameter(1);
% base_param.aR = 0;
% parameter_update = readtable('MCMC_parameter_config_mig1d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'mig1d_1r');
%
% load('../metaData/trait_extraction/gal80d_1r.mat')
% base_param = set_parameter(1);
% base_param.a80 = 0;
% base_param.ag80 = 0;
% parameter_update = readtable('MCMC_parameter_config_gal80d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'gal80d_1r');

load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r.mat')
% parameter_update = readtable('July2nd_parameter_config_wt_set1.csv');

% fit GAL1,3,4; let hill coefficients vary
base_param = set_parameter(1);
parameter_update = readtable('Aug3rd_parameter_config_wt_set1.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'varyN-wildtype_1r', folder1);

% set all hill coefficients equal to 1
base_param = set_parameter(3);
parameter_update = readtable('Aug3rd_parameter_config_wt_set3.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'sequestrate-wildtype_1r', folder2);

%% generate config .mat files that fit one cross

% load('../metaData/trait_extraction/mig1d_1r1c.mat')
% base_param = set_parameter(1);
% base_param.aR = 0;
% parameter_update = readtable('MCMC_parameter_config_mig1d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'mig1d_1r1c');
%
% load('../metaData/trait_extraction/gal80d_1r1c.mat')
% base_param = set_parameter(1);
% base_param.a80 = 0;
% base_param.ag80 = 0;
% parameter_update = readtable('MCMC_parameter_config_gal80d_set1.csv');
% rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'gal80d_1r1c');

load('../metaData/trait_extraction/S288C-double_gradient/wildtype_1r1c.mat')
% parameter_update = readtable('July2nd_parameter_config_wt_set1.csv');

% fit GAL1,3,4; let hill coefficients vary
base_param = set_parameter(1);
parameter_update = readtable('Aug3rd_parameter_config_wt_set1.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'varyN-wildtype_1r1c', folder1);

% set all hill coefficients equal to 1
base_param = set_parameter(3);
parameter_update = readtable('Aug3rd_parameter_config_wt_set3.csv');
rand_init_generator(n_init, trait, n_propose, base_param, parameter_update, 'sequestrate-wildtype_1r1c', folder2);

%% generate shell script for calling mcmc function to run on LSF cluster

% fid = fopen('MCMC_mutant_and_wt_1c_orchestra_Mar30.sh', 'w');
% for i_init = 1:n_init
%     jobtag = 'MCMC_mig1d_1c';
%     filepath = fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']);
%     fprintf(fid, 'bsub -q long -o out -W 240:00 matlab -r "load(''%s''); mcmc2( trait , param_init , parameter_update , ''n_propose'',n_propose , ''error_tol'', error_tol , ''jobtag'', jobtag ); "\n', filepath);
%     jobtag = 'MCMC_gal80d_1c';
%     filepath = fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']);
%     fprintf(fid, 'bsub -q long -o out -W 240:00 matlab -r "load(''%s''); mcmc2( trait , param_init , parameter_update , ''n_propose'',n_propose , ''error_tol'', error_tol , ''jobtag'', jobtag ); "\n', filepath);
%     jobtag = 'MCMC_wt_1c';
%     filepath = fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']);
%     fprintf(fid, 'bsub -q long -o out -W 240:00 matlab -r "load(''%s''); mcmc2( trait , param_init , parameter_update , ''n_propose'',n_propose , ''error_tol'', error_tol , ''jobtag'', jobtag ); "\n', filepath);
% end
% fclose(fid);
