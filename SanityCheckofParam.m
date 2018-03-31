%% part I - mcmc without prior
% load data
mcmc_data_folder = '../results/mcmc_without_prior';
jobtags = {'wildtype_1c', 'mig1d_1c', 'gal80d_1c', ...
    'wildtype_1r', 'mig1d_1r', 'gal80d_1r', ...
    'wildtype_1r1c', 'mig1d_1r1c', 'gal80d_1r1c'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% all update parameters trials
all_update_jobtags = {'all_update_mig1d_1c', 'all_update_gal80d_1c', ...
    'all_update_mig1d_1r', 'all_update_gal80d_1r', ...
    'all_update_mig1d_1r1c', 'all_update_gal80d_1r1c'};
all_update_result = load_mcmc_result(mcmc_data_folder, all_update_jobtags);
% trait = mcmc_result{1,'trait'}{1};
% n_chain = height(mcmc_result);

%% define a threshold for good fitting
load(mcmc_result{1,'filepath'}{1}, 'error_tol');
thresh = - 0.1 / error_tol ^2;  % i.e. obj < 0.1

%% show how parameters fit wildtype data
wt_jobtags = {'wildtype_1c', 'wildtype_1r', 'wildtype_1r1c'};
wt_subset = mcmc_result(ismember(mcmc_result.jobtag, wt_jobtags),:);
wt_subset = sortrows(wt_subset, 'jobtag');
wt_goodfit = wt_subset(wt_subset.param_prob_map > thresh, :);

n_example = height(wt_goodfit);
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
    param = wt_goodfit{i_example, 'param_map'};
    trait = wt_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [wt_conc_Glu_list.row(i_row,:,:), wt_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [wt_conc_Glu_list.col(i_col,:,:), wt_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [wt_conc_Glu_list.cross(i_cross,:,:), wt_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-wildtype-fitting', 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-wildtype-fitting'));

%% show how parameters fit mig1d data
mig1d_jobtags = {'mig1d_1c', 'mig1d_1r', 'mig1d_1r1c'};
mig1d_subset = mcmc_result(ismember(mcmc_result.jobtag, mig1d_jobtags),:);
mig1d_subset = sortrows(mig1d_subset, 'jobtag');
mig1d_goodfit = mig1d_subset(mig1d_subset.param_prob_map > thresh, :);

n_example = height(mig1d_goodfit);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

mig1d_conc_Glu_list = struct();
mig1d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = mig1d_goodfit{i_example, 'param_map'};
    trait = mig1d_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [mig1d_conc_Glu_list.row(i_row,:,:), mig1d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif regexp(mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [mig1d_conc_Glu_list.col(i_col,:,:), mig1d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif regexp(mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [mig1d_conc_Glu_list.cross(i_cross,:,:), mig1d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-mig1d-fitting', 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-mig1d-fitting'));

%% show how parameters fit gal80d data
gal80d_jobtags = {'gal80d_1c', 'gal80d_1r', 'gal80d_1r1c'};
gal80d_subset = mcmc_result(ismember(mcmc_result.jobtag, gal80d_jobtags),:);
gal80d_subset = sortrows(gal80d_subset, 'jobtag');
gal80d_goodfit = gal80d_subset(gal80d_subset.param_prob_map > thresh, :);

n_example = height(gal80d_goodfit);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

gal80d_conc_Glu_list = struct();
gal80d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = gal80d_goodfit{i_example, 'param_map'};
    trait = gal80d_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [gal80d_conc_Glu_list.row(i_row,:,:), gal80d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif regexp(gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [gal80d_conc_Glu_list.col(i_col,:,:), gal80d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif regexp(gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [gal80d_conc_Glu_list.cross(i_cross,:,:), gal80d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-gal80d-fitting', 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-gal80d-fitting'));

%% consolidate all parameters

% all wildtype good fitting parameters
load(wt_subset{1,'filepath'}{1}, 'parameter_update');
param_update_names_wt = parameter_update.parameter_name;
n_update_param_wt = length(param_update_names_wt);
map_param_vals_wt = nan(height(wt_goodfit), n_update_param_wt);

% parameter values in log10 scale
for i_chain = 1:height(wt_goodfit)
    param_map = wt_goodfit{i_chain, 'param_map'};
    for i_param = 1:n_update_param_wt
        map_param_vals_wt(i_chain, i_param) = log10(param_map.(param_update_names_wt{i_param}));
    end
end

n_row = floor(n_update_param_wt ^0.5);
n_col = ceil(n_update_param_wt / n_row);

% all mig1d good fitting parameters
load(mig1d_subset{1,'filepath'}{1}, 'parameter_update');
param_update_names_mig1d = parameter_update.parameter_name;
n_update_param_mig1d = length(param_update_names_mig1d);
map_param_vals_mig1d = nan(height(mig1d_goodfit), n_update_param_mig1d);

% parameter values in log10 scale
for i_chain = 1:height(mig1d_goodfit)
    param_map = mig1d_subset{i_chain, 'param_map'};
    for i_param = 1:n_update_param_mig1d
        map_param_vals_mig1d(i_chain, i_param) = log10(param_map.(param_update_names_mig1d{i_param}));
    end
end

% all gal80d good fitting parameters
load(gal80d_subset{1,'filepath'}{1}, 'parameter_update');
param_update_names_gal80d = parameter_update.parameter_name;
n_update_param_gal80d = length(param_update_names_gal80d);
map_param_vals_gal80d = nan(height(gal80d_goodfit), n_update_param_gal80d);

% parameter values in log10 scale
for i_chain = 1:height(gal80d_goodfit)
    param_map = gal80d_subset{i_chain, 'param_map'};
    for i_param = 1:n_update_param_gal80d
        map_param_vals_gal80d(i_chain, i_param) = log10(param_map.(param_update_names_gal80d{i_param}));
    end
end

%% prepare a color palette
CT = cbrewer('qual', 'Dark2', 8);
c_wt = CT(1,:);
c_mig1d = CT(2,:);
% c_gal80d = CT(3,:);
c_gal80d = 'b';

%% histogram show the difference between wildtype and mig1d
% ksdensity of parameter values of wildtype
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    [f,xi] = ksdensity(map_param_vals_wt(:, i_param));
    plot(xi, f, '-', 'linewidth', 2, 'color', c_wt)
    hold on
%     xlim([min(map_param_vals_wt(:,i_param))-.5, max(map_param_vals_wt(:,i_param))+.5])
    title(param_update_names_wt{i_param})
    grid on
end

% ksdensity of parameter values of mig1d
i_subplot = 1;
for i_param = [1:4,6:n_update_param_wt]     % pass aR
    subplot(n_row, n_col, i_param)
    [f,xi] = ksdensity(map_param_vals_mig1d(:, i_subplot));
    plot(xi, f, '-', 'linewidth', 2, 'color', c_mig1d)
    hold on
    i_subplot = i_subplot+1;
end

% show parameter prior distribution
load(wt_subset{1,'filepath'}{1}, 'parameter_update');
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    prior_mean = parameter_update{i_param, 'prior_mean'};
    prior_sigma = parameter_update{i_param, 'prior_sigma'};
    xi = xlim();
    xi = linspace(xi(1), xi(2), 1000);
    yi = pdf('normal', xi, log10(prior_mean), prior_sigma * log10(exp(1)));
    yi = yi/max(yi) * max(ylim()) * 0.8;  % scaling so that easier to compare
    plot(xi, yi, 'k:', 'linewidth', 2)
end

suplabel(changeunderscore('wt_compare_with_mig1d_filtered'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'wt_compare_with_mig1d_filtered'));

%% histogram show the difference between wildtype and gal80d
% ksdensity of parameter values of wildtype
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    [f,xi] = ksdensity(map_param_vals_wt(:, i_param));
    plot(xi, f, '-', 'linewidth', 2, 'color', c_wt)
    hold on
%     xlim([min(map_param_vals_wt(:,i_param))-.5, max(map_param_vals_wt(:,i_param))+.5])
    title(param_update_names_wt{i_param})
    grid on
end

% ksdensity of parameter values of gal80d
i_subplot = 1;
for i_param = [1:3, 5:8, 10:n_update_param_wt]     % pass a80, ag80
    subplot(n_row, n_col, i_param)
    [f,xi] = ksdensity(map_param_vals_gal80d(:, i_subplot));
    plot(xi, f, '-', 'linewidth', 2, 'color', c_gal80d)
    hold on
    i_subplot = i_subplot+1;
end

% show parameter prior distribution
load(wt_subset{1,'filepath'}{1}, 'parameter_update');
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    prior_mean = parameter_update{i_param, 'prior_mean'};
    prior_sigma = parameter_update{i_param, 'prior_sigma'};
    xi = xlim();
    xi = linspace(xi(1), xi(2), 1000);
    yi = pdf('normal', xi, log10(prior_mean), prior_sigma * log10(exp(1)));
    yi = yi/max(yi) * max(ylim()) * 0.8;  % scaling so that easier to compare
    plot(xi, yi, 'k:', 'linewidth', 2)
end

suplabel(changeunderscore('wt_compare_with_gal80d_filtered'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'wt_compare_with_gal80d_filtered'));

%% scatter plot of the parameter values across wt, mig1d and gal80d
xlimub = max([size(map_param_vals_gal80d,1), size(map_param_vals_mig1d,1), size(map_param_vals_wt,1)]);
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_wt(:, i_param), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 xlimub+1])
    title(param_update_names_wt{i_param})
    grid on
end

i_subplot = 1;
for i_param = [1:4,6:n_update_param_wt]     % pass aR
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_mig1d(:, i_subplot), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
    i_subplot = i_subplot+1;
end

i_subplot = 1;
for i_param = [1:3, 5:8, 10:n_update_param_wt]     % pass a80, ag80
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_gal80d(:, i_subplot), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
    i_subplot = i_subplot+1;
end

if i_param == n_update_param_wt
    h = legend('wt', 'mig1d', 'gal80d', 'location', 'e');
    tmp = get(h, 'position');
    tmp(1) = tmp(1) + .2;
    set(h, 'position', tmp);
end

suplabel(changeunderscore('scatter_of_wt_and_mutant_filtered'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'scatter_of_wt_and_mutant_filtered'));

%% show the steady state concentrations for 12 species
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1', 'Mig1*', 'C83', 'C84', 'glu', 'gal'};
n_species = length(species_list);
n_row = floor(n_species ^0.5);
n_col = ceil(n_species / n_row);

% plot the steady state concentrations for one column, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], wt_conc_Gal_list.col(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 9])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], mig1d_conc_Gal_list.col(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], gal80d_conc_Gal_list.col(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

suplabel(changeunderscore('S.S.conc_one_column'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_column'));

%% plot the steady state concentrations for one row, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], wt_conc_Gal_list.row(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 13])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], mig1d_conc_Gal_list.row(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], gal80d_conc_Gal_list.row(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

suplabel(changeunderscore('S.S.conc_one_row'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_row'));

% plot the steady state concentrations for one cross, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:19], wt_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 20])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:19], mig1d_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:19], gal80d_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

suplabel(changeunderscore('S.S.conc_one_cross'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_cross'));

%% sanity check for all parameters updated in mutants
% i.e. not forcing aR=0 in mig1d && not forcing a80=ag80=0 in gal80d

% show how parameters fit mig1d data without forcing aR=0
all_update_mig1d_jobtags = {'all_update_mig1d_1c', 'all_update_mig1d_1r', 'all_update_mig1d_1r1c'};
all_update_mig1d_subset = all_update_result(ismember(all_update_result.jobtag, all_update_mig1d_jobtags),:);
all_update_mig1d_subset = sortrows(all_update_mig1d_subset, 'jobtag');
all_update_mig1d_goodfit = all_update_mig1d_subset(all_update_mig1d_subset.param_prob_map > thresh, :);

n_example = height(all_update_mig1d_goodfit);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

all_update_mig1d_conc_Glu_list = struct();
all_update_mig1d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = all_update_mig1d_goodfit{i_example, 'param_map'};
    trait = all_update_mig1d_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(all_update_mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [all_update_mig1d_conc_Glu_list.row(i_row,:,:), all_update_mig1d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif regexp(all_update_mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [all_update_mig1d_conc_Glu_list.col(i_col,:,:), all_update_mig1d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif regexp(all_update_mig1d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [all_update_mig1d_conc_Glu_list.cross(i_cross,:,:), all_update_mig1d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel(changeunderscore('all_update_filtered-mig1d-fitting'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'all_update_filtered-mig1d-fitting'));

% show how parameters fit gal80d data without forcing a80=ag80=0
all_update_gal80d_jobtags = {'all_update_gal80d_1c', 'all_update_gal80d_1r', 'all_update_gal80d_1r1c'};
all_update_gal80d_subset = all_update_result(ismember(all_update_result.jobtag, all_update_gal80d_jobtags),:);
all_update_gal80d_subset = sortrows(all_update_gal80d_subset, 'jobtag');
all_update_gal80d_goodfit = all_update_gal80d_subset(all_update_gal80d_subset.param_prob_map > thresh, :);

n_example = height(all_update_gal80d_goodfit);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

all_update_gal80d_conc_Glu_list = struct();
all_update_gal80d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = all_update_gal80d_goodfit{i_example, 'param_map'};
    trait = all_update_gal80d_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)

    % mapping fit type with jobtags
    if regexp(all_update_gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [all_update_gal80d_conc_Glu_list.row(i_row,:,:), all_update_gal80d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif regexp(all_update_gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [all_update_gal80d_conc_Glu_list.col(i_col,:,:), all_update_gal80d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif regexp(all_update_gal80d_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [all_update_gal80d_conc_Glu_list.cross(i_cross,:,:), all_update_gal80d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel(changeunderscore('all_update_filtered-gal80d-fitting'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'all_update_filtered-gal80d-fitting'));

% consolidate all parameters

% all mig1d good fitting parameters without forcing aR=0
load(all_update_mig1d_subset{1,'filepath'}{1}, 'parameter_update');
param_update_names_mig1d = parameter_update.parameter_name;
n_update_param_mig1d = length(param_update_names_mig1d);
map_param_vals_mig1d_all_update = nan(height(all_update_mig1d_goodfit), n_update_param_mig1d);

% parameter values in log10 scale
for i_chain = 1:height(all_update_mig1d_goodfit)
    param_map = all_update_mig1d_subset{i_chain, 'param_map'};
    for i_param = 1:n_update_param_mig1d
        map_param_vals_mig1d_all_update(i_chain, i_param) = log10(param_map.(param_update_names_mig1d{i_param}));
    end
end

% all gal80d good fitting parameters without forcing a80=ag80=0
load(all_update_gal80d_subset{1,'filepath'}{1}, 'parameter_update');
param_update_names_gal80d = parameter_update.parameter_name;
n_update_param_gal80d = length(param_update_names_gal80d);
map_param_vals_gal80d_all_update = nan(height(all_update_gal80d_goodfit), n_update_param_gal80d);

% parameter values in log10 scale
for i_chain = 1:height(all_update_gal80d_goodfit)
    param_map = all_update_gal80d_subset{i_chain, 'param_map'};
    for i_param = 1:n_update_param_gal80d
        map_param_vals_gal80d_all_update(i_chain, i_param) = log10(param_map.(param_update_names_gal80d{i_param}));
    end
end

% scatter plot of the parameter values across wt, mig1d and gal80d - all parameters update in mutants
n_row = floor(n_update_param_wt ^0.5);
n_col = ceil(n_update_param_wt / n_row);
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_wt(:, i_param), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 size(map_param_vals_wt,1)])
    title(param_update_names_wt{i_param})
    grid on
end

for i_param = 1:n_update_param_wt    
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_mig1d_all_update(:, i_param), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end

for i_param = 1:n_update_param_wt  
    subplot(n_row, n_col, i_param)
    plot(map_param_vals_gal80d_all_update(:, i_param), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

if i_param == n_update_param_wt
    h = legend('wt', 'mig1d', 'gal80d', 'location', 'e');
    tmp = get(h, 'position');
    tmp(1) = tmp(1) + .2;
    set(h, 'position', tmp);
end

suplabel(changeunderscore('scatter_of_wt_and_mutant_filtered-all-update'), 't');
export_fig(fullfile('../results/param_sanity_check_plot/', 'scatter_of_wt_and_mutant_filtered-all-update'));

% show the steady state concentrations for 12 species - all parameters update in mutants
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1', 'Mig1*', 'C83', 'C84', 'glu', 'gal'};
n_species = length(species_list);
n_row = floor(n_species ^0.5);
n_col = ceil(n_species / n_row);

% plot the steady state concentrations for one column, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], wt_conc_Gal_list.col(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 9])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], all_update_mig1d_conc_Gal_list.col(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:8], all_update_gal80d_conc_Gal_list.col(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

suplabel(changeunderscore('S.S.conc_one_column_all_update'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_column_all_update'));

% plot the steady state concentrations for one row, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], wt_conc_Gal_list.row(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 13])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], all_update_mig1d_conc_Gal_list.row(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:12], all_update_gal80d_conc_Gal_list.row(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
    hold on
end

suplabel(changeunderscore('S.S.conc_one_row_all_update'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_row_all_update'));

% plot the steady state concentrations for one cross, only plot ON state
figure
set(gcf, 'position', [298 107 1259 822])
% plot wild-type
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:19], wt_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_wt, 'markersize', 15)
    hold on
    xlim([0 20])
    set(gca, 'yscale', 'log')
    title(species_list{i_species})
    grid on
end
% plot mig1d
for i_species = 1:n_species
    subplot(n_row, n_col, i_species)
    plot([1:19], all_update_mig1d_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_mig1d, 'markersize', 15)
    hold on
end
% plot gal80d
if isfield(all_update_gal80d_conc_Gal_list, 'cross')
    for i_species = 1:n_species
        subplot(n_row, n_col, i_species)
        plot([1:19], all_update_gal80d_conc_Gal_list.cross(:,:,i_species), '.', 'color', c_gal80d, 'markersize', 15)
        hold on
    end
end

suplabel(changeunderscore('S.S.conc_one_cross_all_update'), 't');
% export figure
export_fig(fullfile('../results/param_sanity_check_plot/', 'S.S.conc_one_cross_all_update'));


