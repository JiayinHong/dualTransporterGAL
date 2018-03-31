%% part I - sanity check of the parameters optimized by 'simulateGALPathway'
% load data
mcmc_data_folder = '../results/mcmc_multiple_trait';
jobtags = {'triple_fit_1c', 'triple_fit_1r', 'triple_fit_1r1c'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% define a threshold for good fitting
% load(mcmc_result{1,'filepath'}{1}, 'error_tol');
% thresh = - 0.1 / error_tol ^2;  % i.e. obj < 0.1

%% show how parameters fit wildtype data
n_example = height(mcmc_result);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

wt_conc_Glu_list = struct();
wt_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = mcmc_result{i_example, 'param_map'};
    trait = mcmc_result{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [wt_conc_Glu_list.row(i_row,:,:), wt_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait.wt, fit_type);
        i_row = i_row + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [wt_conc_Glu_list.col(i_col,:,:), wt_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait.wt, fit_type);
        i_col = i_col + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [wt_conc_Glu_list.cross(i_cross,:,:), wt_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait.wt, fit_type);
        i_cross = i_cross + 1;
    end
    
end
% suplabel('filtered-wildtype-fitting', 't');
% export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-wildtype-fitting'));

%% show how parameters fit mig1d data
mig1d_conc_Glu_list = struct();
mig1d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = mcmc_result{i_example, 'param_map'};
    param.aR = 0;
    trait = mcmc_result{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [mig1d_conc_Glu_list.row(i_row,:,:), mig1d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait.mig1d, fit_type);
        i_row = i_row + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [mig1d_conc_Glu_list.col(i_col,:,:), mig1d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait.mig1d, fit_type);
        i_col = i_col + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [mig1d_conc_Glu_list.cross(i_cross,:,:), mig1d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait.mig1d, fit_type);
        i_cross = i_cross + 1;
    end
    
end
% suplabel('filtered-mig1d-fitting', 't');
% export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-mig1d-fitting'));

%% show how parameters fit gal80d data
gal80d_conc_Glu_list = struct();
gal80d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = mcmc_result{i_example, 'param_map'};
    param.a80 = 0;
    param.ag80 = 0;
    trait = mcmc_result{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [gal80d_conc_Glu_list.row(i_row,:,:), gal80d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait.gal80d, fit_type);
        i_row = i_row + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [gal80d_conc_Glu_list.col(i_col,:,:), gal80d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait.gal80d, fit_type);
        i_col = i_col + 1;
    elseif regexp(mcmc_result{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [gal80d_conc_Glu_list.cross(i_cross,:,:), gal80d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait.gal80d, fit_type);
        i_cross = i_cross + 1;
    end
    
end
% suplabel('filtered-gal80d-fitting', 't');
% export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-gal80d-fitting'));


