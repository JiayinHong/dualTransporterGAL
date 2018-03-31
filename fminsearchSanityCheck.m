clear
clc

fminsearch_data_folder = '../results/fminsearch/';
files = dir( fullfile(fminsearch_data_folder, '*.txt') );
n_file = length(files);
filename_pat = '([\w_\d]+)-\d{3}-\d{6}_\d{2}:\d{2}.txt';
% example: 'wildtype_1r1c-003-170502_22:10.txt'
jobtags = {'wildtype_1c', 'mig1d_1c', 'gal80d_1c', ...
    'wildtype_1r', 'mig1d_1r', 'gal80d_1r', ...
    'wildtype_1r1c', 'mig1d_1r1c', 'gal80d_1r1c'};

thresh_obj = 0.1;

i_wt = 1;
i_mig1d = 1;
i_gal80d = 1;
for i_file = 1:n_file
    filename = files(i_file).name;
    tok = regexp(filename, filename_pat, 'tokens');
    if isempty(tok)
        fprintf('File name format warning: %s\n', filename)
    elseif regexp(tok{1}{1}, 'wildtype')
        filepath = fullfile(fminsearch_data_folder, filename);
        wt_array{i_wt} = filepath;
        i_wt = i_wt + 1;
    elseif regexp(tok{1}{1}, 'mig1d')
        filepath = fullfile(fminsearch_data_folder, filename);
        mig1d_array{i_mig1d} = filepath;
        i_mig1d = i_mig1d + 1;
    elseif regexp(tok{1}{1}, 'gal80d')
        filepath = fullfile(fminsearch_data_folder, filename);
        gal80d_array{i_gal80d} = filepath;
        i_gal80d = i_gal80d + 1;
    end
end

wt_array = wt_array';
mig1d_array = mig1d_array';
gal80d_array = gal80d_array';

%% show how parameters fit wildtype data
base_param = set_parameter(1);
wt_param_update = readtable('MCMC_parameter_config_wt_set1.csv');
parameter_name = wt_param_update.parameter_name;
load('../metaData/trait_extraction/wildtype_1r.mat');
wt_trait_1r = trait;
load('../metaData/trait_extraction/wildtype_1c.mat');
wt_trait_1c = trait;
load('../metaData/trait_extraction/wildtype_1r1c.mat');
wt_trait_1r1c = trait;

n_wt = length(wt_array);
param_list = cell(n_wt,1);
obj_list = nan(n_wt,1);
iter_list = nan(n_wt,1);
fit_type_list = cell(n_wt,1);
trait_list = cell(n_wt,1);
for i = 1:n_wt
    [param_values, obj, iter] = read_fminsearch_result(wt_array{i});
    param_list{i} = update_param(base_param, parameter_name, param_values);
    obj_list(i) = obj;
    iter_list(i) = iter;
    if regexp(wt_array{i}, '.*_1r-.*')
        fit_type_list{i} = 'one_row';
        trait_list{i} = wt_trait_1r;
    elseif regexp(wt_array{i}, '.*_1c-.*')
        fit_type_list{i} = 'one_column';
        trait_list{i} = wt_trait_1c;
    elseif regexp(wt_array{i}, '.*_1r1c-.*')
        fit_type_list{i} = 'one_cross';
        trait_list{i} = wt_trait_1r1c;
    end
end
wt_fmin_tab = table(wt_array, iter_list, obj_list, ...
    param_list, fit_type_list, trait_list, ...
    'VariableNames', {'filepath', 'iteration', 'obj', ...
    'param', 'fit_type', 'trait'});
wt_goodfit = wt_fmin_tab(wt_fmin_tab.obj < thresh_obj, :);

n_example = height(wt_goodfit);
n_row = floor(n_example ^.5);
n_col = ceil(n_example / n_row);

wt_conc_Glu_list = struct();
wt_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 622 1217 290]);
for i_example = 1:n_example
    param = wt_goodfit{i_example, 'param'}{1};
    trait = wt_goodfit{i_example, 'trait'}{1};
    fit_type = wt_goodfit{i_example, 'fit_type'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if strcmp(fit_type, 'one_row')
        [wt_conc_Glu_list.row(i_row,:,:), wt_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif strcmp(fit_type, 'one_column')
        [wt_conc_Glu_list.col(i_col,:,:), wt_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif strcmp(fit_type, 'one_cross')
        [wt_conc_Glu_list.cross(i_cross,:,:), wt_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-wildtype-fitting', 't');
export_fig(fullfile('../results/fminsearch_sanity_check_plot/', 'filtered-wildtype-fitting'));

%% show how parameters fit mig1d data
base_param = set_parameter(1);
base_param.aR = 0;
mig1d_param_update = readtable('MCMC_parameter_config_mig1d_set1.csv');
parameter_name = mig1d_param_update.parameter_name;
load('../metaData/trait_extraction/mig1d_1r.mat');
mig1d_trait_1r = trait;
load('../metaData/trait_extraction/mig1d_1c.mat');
mig1d_trait_1c = trait;
load('../metaData/trait_extraction/mig1d_1r1c.mat');
mig1d_trait_1r1c = trait;

n_mig1d = length(mig1d_array);
param_list = cell(n_mig1d,1);
obj_list = nan(n_mig1d,1);
iter_list = nan(n_mig1d,1);
fit_type_list = cell(n_mig1d,1);
trait_list = cell(n_mig1d,1);
for i = 1:n_mig1d
    [param_values, obj, iter] = read_fminsearch_result(mig1d_array{i});
    param_list{i} = update_param(base_param, parameter_name, param_values);
    obj_list(i) = obj;
    iter_list(i) = iter;
    if regexp(mig1d_array{i}, '.*_1r-.*')
        fit_type_list{i} = 'one_row';
        trait_list{i} = mig1d_trait_1r;
    elseif regexp(mig1d_array{i}, '.*_1c-.*')
        fit_type_list{i} = 'one_column';
        trait_list{i} = mig1d_trait_1c;
    elseif regexp(mig1d_array{i}, '.*_1r1c-.*')
        fit_type_list{i} = 'one_cross';
        trait_list{i} = mig1d_trait_1r1c;
    end
end
mig1d_fmin_tab = table(mig1d_array, iter_list, obj_list, ...
    param_list, fit_type_list, trait_list, ...
    'VariableNames', {'filepath', 'iteration', 'obj', ...
    'param', 'fit_type', 'trait'});
mig1d_goodfit = mig1d_fmin_tab(mig1d_fmin_tab.obj < thresh_obj, :);

n_example = height(mig1d_goodfit);
n_row = floor(n_example ^.5);
n_col = ceil(n_example / n_row);

mig1d_conc_Glu_list = struct();
mig1d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 622 1217 290]);
for i_example = 1:n_example
    param = mig1d_goodfit{i_example, 'param'}{1};
    trait = mig1d_goodfit{i_example, 'trait'}{1};
    fit_type = mig1d_goodfit{i_example, 'fit_type'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if strcmp(fit_type, 'one_row')
        [mig1d_conc_Glu_list.row(i_row,:,:), mig1d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif strcmp(fit_type, 'one_column')
        [mig1d_conc_Glu_list.col(i_col,:,:), mig1d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif strcmp(fit_type, 'one_cross')
        [mig1d_conc_Glu_list.cross(i_cross,:,:), mig1d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-mig1d-fitting', 't');
export_fig(fullfile('../results/fminsearch_sanity_check_plot/', 'filtered-mig1d-fitting'));

%% show how parameters fit gal80d data
base_param = set_parameter(1);
base_param.a80 = 0;
base_param.ag80 = 0;
gal80d_param_update = readtable('MCMC_parameter_config_gal80d_set1.csv');
parameter_name = gal80d_param_update.parameter_name;
load('../metaData/trait_extraction/gal80d_1r.mat');
gal80d_trait_1r = trait;
load('../metaData/trait_extraction/gal80d_1c.mat');
gal80d_trait_1c = trait;
load('../metaData/trait_extraction/gal80d_1r1c.mat');
gal80d_trait_1r1c = trait;

n_gal80d = length(gal80d_array);
param_list = cell(n_gal80d,1);
obj_list = nan(n_gal80d,1);
iter_list = nan(n_gal80d,1);
fit_type_list = cell(n_gal80d,1);
trait_list = cell(n_gal80d,1);
for i = 1:n_gal80d
    [param_values, obj, iter] = read_fminsearch_result(gal80d_array{i});
    param_list{i} = update_param(base_param, parameter_name, param_values);
    obj_list(i) = obj;
    iter_list(i) = iter;
    if regexp(gal80d_array{i}, '.*_1r-.*')
        fit_type_list{i} = 'one_row';
        trait_list{i} = gal80d_trait_1r;
    elseif regexp(gal80d_array{i}, '.*_1c-.*')
        fit_type_list{i} = 'one_column';
        trait_list{i} = gal80d_trait_1c;
    elseif regexp(gal80d_array{i}, '.*_1r1c-.*')
        fit_type_list{i} = 'one_cross';
        trait_list{i} = gal80d_trait_1r1c;
    end
end
gal80d_fmin_tab = table(gal80d_array, iter_list, obj_list, ...
    param_list, fit_type_list, trait_list, ...
    'VariableNames', {'filepath', 'iteration', 'obj', ...
    'param', 'fit_type', 'trait'});
gal80d_goodfit = gal80d_fmin_tab(gal80d_fmin_tab.obj < thresh_obj, :);

n_example = height(gal80d_goodfit);
n_row = floor(n_example ^.5);
n_col = ceil(n_example / n_row);

gal80d_conc_Glu_list = struct();
gal80d_conc_Gal_list = struct();
i_row = 1;
i_col = 1;
i_cross = 1;

figure
set(gcf, 'position', [401 407 1264 505]);
for i_example = 1:n_example
    param = gal80d_goodfit{i_example, 'param'}{1};
    trait = gal80d_goodfit{i_example, 'trait'}{1};
    fit_type = gal80d_goodfit{i_example, 'fit_type'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if strcmp(fit_type, 'one_row')
        [gal80d_conc_Glu_list.row(i_row,:,:), gal80d_conc_Gal_list.row(i_row,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_row = i_row + 1;
    elseif strcmp(fit_type, 'one_column')
        [gal80d_conc_Glu_list.col(i_col,:,:), gal80d_conc_Gal_list.col(i_col,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_col = i_col + 1;
    elseif strcmp(fit_type, 'one_cross')
        [gal80d_conc_Glu_list.cross(i_cross,:,:), gal80d_conc_Gal_list.cross(i_cross,:,:)] = param_fitting_plot(param, trait, fit_type);
        i_cross = i_cross + 1;
    end
    
end
suplabel('filtered-gal80d-fitting', 't');
export_fig(fullfile('../results/fminsearch_sanity_check_plot/', 'filtered-gal80d-fitting'));
