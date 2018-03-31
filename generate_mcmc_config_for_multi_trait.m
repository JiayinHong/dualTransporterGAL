%%
n_init = 200;   % how many random initial states
n_propose = 50000;     % how many iterations
error_tol = .15;

folder_random_init = '../metaData/mutant_and_wt_triple_fit';
if ~exist(folder_random_init)
    mkdir(folder_random_init)
end

base_param = set_parameter(8);
parameter_update = readtable('Mar3rd_param_config_set8.csv');

%% generate wt/mig1d/gal80d triple fit config files 
% updated by JH on March 4th, 2018

trait_multiple = struct();
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat')
trait_multiple.wildtype = trait;
load('../metaData/trait_extraction/S288C-double_gradient/mig1d_all_data.mat')
trait_multiple.mig1d = trait;
load('../metaData/trait_extraction/S288C-double_gradient/gal80d_all_data.mat')
trait_multiple.gal80d = trait;

for i_init = 601:n_init+600
    parameter_val = nan(height(parameter_update),1);
    for i = 1:height(parameter_update)
        parameter_val(i) = random('lognormal', ...
            log(parameter_update{i, 'prior_mean'}), ...
            parameter_update{i, 'prior_sigma'} ...
            );
    end
    param_init = update_param(base_param, parameter_update.parameter_name, parameter_val);
    
    trait = trait_multiple;
    jobtag = 'triple_fit';
    save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']) ...
        , 'parameter_update', 'param_init' ...
        , 'error_tol', 'n_propose' ...
        , 'trait' ...
        , 'jobtag' ...
        )
end

% parameter_update_wt = readtable('MCMC_parameter_config_wt.csv');
% parameter_update_mig1d = readtable('MCMC_parameter_config_mig1d.csv');
% parameter_update_gal80d = readtable('MCMC_parameter_config_gal80d.csv');
% parameter_update = struct();
% parameter_update.wt = parameter_update_wt;
% parameter_update.mig1d = parameter_update_mig1d;
% parameter_update.gal80d = parameter_update_gal80d;
% param_init_base_wt = set_parameter(2);

%% generate wt/mig1d/gal80d triple fit config files for one column
trait_multiple = struct();
load('../metaData/trait_extraction/wildtype_1c.mat')
trait_multiple.wt = trait;
load('../metaData/trait_extraction/mig1d_1c.mat')
trait_multiple.mig1d = trait;
load('../metaData/trait_extraction/gal80d_1c.mat')
trait_multiple.gal80d = trait;

for i_init = 1:n_init
    parameter_val_wt = nan(height(parameter_update_wt),1);
    for i = 1:height(parameter_update_wt)
        parameter_val_wt(i) = random('lognormal', ...
            log(parameter_update_wt{i, 'prior_mean'}), ...
            parameter_update_wt{i, 'prior_sigma'} ...
            );
    end
    param_init_wt = update_param(param_init_base_wt, parameter_update_wt.parameter_name, parameter_val_wt);
    
    param_init_mig1d = param_init_wt;
    param_init_mig1d.aR = 0;
    
    param_init_gal80d = param_init_wt;
    param_init_gal80d.a80 = 0;
    param_init_gal80d.ag80 = 0;
    
    param_init = struct();
    param_init.wt = param_init_wt;
    param_init.mig1d = param_init_mig1d;
    param_init.gal80d = param_init_gal80d;
    
    trait = trait_multiple;
    jobtag = 'triple_fit_1c';
    save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']) ...
        , 'parameter_update', 'param_init' ...
        , 'error_tol', 'n_propose' ...
        , 'trait' ...
        , 'jobtag' ...
        )
end

%% generate wt/mig1d/gal80d triple fit config files for one row
trait_multiple = struct();
load('../metaData/trait_extraction/wildtype_1r.mat')
trait_multiple.wt = trait;
load('../metaData/trait_extraction/mig1d_1r.mat')
trait_multiple.mig1d = trait;
load('../metaData/trait_extraction/gal80d_1r.mat')
trait_multiple.gal80d = trait;

for i_init = 1:n_init
    parameter_val_wt = nan(height(parameter_update_wt),1);
    for i = 1:height(parameter_update_wt)
        parameter_val_wt(i) = random('lognormal', ...
            log(parameter_update_wt{i, 'prior_mean'}), ...
            parameter_update_wt{i, 'prior_sigma'} ...
            );
    end
    param_init_wt = update_param(param_init_base_wt, parameter_update_wt.parameter_name, parameter_val_wt);
    
    param_init_mig1d = param_init_wt;
    param_init_mig1d.aR = 0;
    
    param_init_gal80d = param_init_wt;
    param_init_gal80d.a80 = 0;
    param_init_gal80d.ag80 = 0;
    
    param_init = struct();
    param_init.wt = param_init_wt;
    param_init.mig1d = param_init_mig1d;
    param_init.gal80d = param_init_gal80d;
    
    trait = trait_multiple;
    jobtag = 'triple_fit_1r';
    save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']) ...
        , 'parameter_update', 'param_init' ...
        , 'error_tol', 'n_propose' ...
        , 'trait' ...
        , 'jobtag' ...
        )
end

%% generate wt/mig1d/gal80d triple fit config files for one cross
trait_multiple = struct();
load('../metaData/trait_extraction/wildtype_1r1c.mat')
trait_multiple.wt = trait;
load('../metaData/trait_extraction/mig1d_1r1c.mat')
trait_multiple.mig1d = trait;
load('../metaData/trait_extraction/gal80d_1r1c.mat')
trait_multiple.gal80d = trait;

for i_init = 1:n_init
    parameter_val_wt = nan(height(parameter_update_wt),1);
    for i = 1:height(parameter_update_wt)
        parameter_val_wt(i) = random('lognormal', ...
            log(parameter_update_wt{i, 'prior_mean'}), ...
            parameter_update_wt{i, 'prior_sigma'} ...
            );
    end
    param_init_wt = update_param(param_init_base_wt, parameter_update_wt.parameter_name, parameter_val_wt);
    
    param_init_mig1d = param_init_wt;
    param_init_mig1d.aR = 0;
    
    param_init_gal80d = param_init_wt;
    param_init_gal80d.a80 = 0;
    param_init_gal80d.ag80 = 0;
    
    param_init = struct();
    param_init.wt = param_init_wt;
    param_init.mig1d = param_init_mig1d;
    param_init.gal80d = param_init_gal80d;
    
    trait = trait_multiple;
    jobtag = 'triple_fit_1r1c';
    save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat']) ...
        , 'parameter_update', 'param_init' ...
        , 'error_tol', 'n_propose' ...
        , 'trait' ...
        , 'jobtag' ...
        )
end
