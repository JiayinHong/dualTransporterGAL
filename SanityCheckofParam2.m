% this script is modified from 'SanityCheckofParam'
% 2017.07.10
% the original script 'SanityCheckofParam' includes the basic code for
% making plots, was created to analyze 'wildtype', 'mig1d', and 'gal80d'
% strains, with forcing knockout corresponding parameters = 0, or with all
% parameters to be varied, called 'all update'.

% this second version is modified to dealing with mcmc results that using
% GAL3pr-YFP, GAL4pr-YFP data to further constrain the model. Applied to
% wildtype data only.

%% fit GAL1, 3 & 4 - piror included - external glucose replace R*
mcmc_data_folder = '../results/external_glucose_replace_R*/';
jobtags = {'medium-wildtype_1c', 'medium-wildtype_1r'};
mcmc_old_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, 3 & 4 - prior included - remove R*
mcmc_data_folder = '../results/mcmc_for_GAL234-removeR/';
jobtags = {'medium-wildtype_1r', 'medium-wildtype_1r1c'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, 3 & 4 - prior included - change the form of R*
% medium step size
mcmc_data_folder = '../results/mutants&GAL34/';
jobtags = {'medium-wildtype_1r1c', 'medium-gal80d_1c', 'medium-gal80d_1r'...
          , 'medium-mig1d_1c', 'medium-mig1d_1r'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, 3 & 4 - prior included - change the form of R*
% medium step size
mcmc_data_folder = '../results/changedRform/';
jobtags = {'medium-wildtype_1r', 'medium-wildtype_1c'};
mcmc_result_medium_step = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, 3 & 4 - prior included - use alpha*KMglu to replace KMgal
% also use beta*kglu to replace kgal
% medium step size
mcmc_data_folder = '../results/fitGAL134-MediumStepSize/';
jobtags = {'medium-wildtype_1r', 'medium-wildtype_1c', 'medium-wildtype_1r1c'};
mcmc_result_medium_step = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, 3 & 4 - prior included - use alpha*KMglu to replace KMgal
% also use beta*kglu to replace kgal
% test step size and span of prior distribution
mcmc_data_folder = '../results/mcmc_for_GAL234-test_stepsize/';
jobtags = {'small-wildtype_1r', 'small-wildtype_1c', 'small-wildtype_1r1c'...
              ,'medium-wildtype_1r', 'medium-wildtype_1c', 'medium-wildtype_1r1c'...
              ,'large-wildtype_1r', 'large-wildtype_1c', 'large-wildtype_1r1c'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, GAL3 & GAL4 - prior included version
% large step size
mcmc_data_folder = '../rvesults/LargeStepSize-mcmc_for_GAL234-prior_included';
jobtags = {'varyN-wildtype_1c', 'varyN-wildtype_1r', 'varyN-wildtype_1r1c'...
            , 'sequestrate-wildtype_1c', 'sequestrate-wildtype_1r', 'sequestrate-wildtype_1r1c'};
mcmc_result_large_step = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, GAL3 & GAL4 - prior included version
% medium step size
mcmc_data_folder = '../results/MediumStepSize-mcmc_for_GAL234-prior_included';
jobtags = {'varyN-wildtype_1c', 'varyN-wildtype_1r', 'varyN-wildtype_1r1c'...
            , 'sequestrate-wildtype_1c', 'sequestrate-wildtype_1r', 'sequestrate-wildtype_1r1c'};
mcmc_result_medium_step = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, GAL3 & GAL4 - prior included version
% small step size
mcmc_data_folder = '../results/SmallStepSize-mcmc_for_GAL234-prior_included';
jobtags = {'wildtype_1c', 'wildtype_1r', 'wildtype_1r1c'};
mcmc_result_small_step = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, GAL3 & GAL4
mcmc_data_folder = '../results/mcmc_for_GAL234';
jobtags = {'wildtype_1c', 'wildtype_1r', 'wildtype_1r1c'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% fit GAL1, GAL3 & GAL4, long chains - 1000,000 iterations
mcmc_data_folder = '../results/GAL234-longChains/';
jobtags = {'wildtype_1c', 'wildtype_1r', 'wildtype_1r1c'...
            , 'MAP-wildtype_1c', 'MAP-wildtype_1r', 'MAP-wildtype_1r1c'...
            };
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% define a threshold for good fitting
load(mcmc_result{1,'filepath'}{1}, 'error_tol');
thresh = - 0.1 / error_tol ^2;  % i.e. obj < 0.1
% set font size and marker size
markersize = 10;
fontsize = 12;

%% first, show how parameters fit wildtype GAL1 level
sorted_wt = sortrows(mcmc_result, 'jobtag');
wt_goodfit = sorted_wt(sorted_wt.param_prob_map > thresh, :);

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
[~,h] = suplabel('filtered-wildtype-fitting', 't');
h.FontSize = fontsize;
% export_fig(fullfile('../results/param_sanity_check_plot/', 'filtered-wildtype-fitting'));

%% Then, show how the parameters fit GAL3, GAL4 level
load('../metaData/trait_extraction/GAL3pr_all_data.mat')
G3level = trait.basal_level;
load('../metaData/trait_extraction/GAL4pr_all_data.mat')
G4level = trait.basal_level;

% plot G3
FLAG = 'G3';
figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = wt_goodfit{i_example, 'param_map'};
    trait = wt_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    end
    
end
[~,h] = suplabel('GAL3 level', 't');
h.FontSize = fontsize;

% plot G4
FLAG = 'G4';
figure
set(gcf, 'position', [401 155 1262 757]);
for i_example = 1:n_example
    param = wt_goodfit{i_example, 'param_map'};
    trait = wt_goodfit{i_example, 'trait'}{1};
    subplot(n_row, n_col, i_example)
    
    % mapping fit type with jobtags
    if regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    elseif regexp(wt_goodfit{i_example, 'jobtag'}{1}, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
        [~,~] = param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, FLAG);
    end  
    
end
[~,h] = suplabel('GAL4 level', 't');
h.FontSize = fontsize;

%% all wildtype good fitting parameters
load(sorted_wt{1,'filepath'}{1}, 'parameter_update');
param_update_names_wt = parameter_update.parameter_name;
n_update_param_wt = length(param_update_names_wt);
map_param_vals_wt = nan(height(wt_goodfit), n_update_param_wt);

% parameter values in log10 scale
for i_chain = 1:height(wt_goodfit)
    param_map = wt_goodfit{i_chain, 'param_map'};
    for i_param = 1:n_update_param_wt
        map_param_vals_wt(i_chain, i_param) = param_map.(param_update_names_wt{i_param});
    end
end

n_row = floor(n_update_param_wt ^0.5);
n_col = ceil(n_update_param_wt / n_row);

%% scatter of parameter values of wildtype
n_chain = size(map_param_vals_wt,1);
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    for i_chain = 1:n_chain
        plot(i_chain, map_param_vals_wt(i_chain, i_param), '.', 'markersize', markersize);
        hold on
    end
    
    xlim([0 n_chain+1])
    set(gca, 'yscale', 'log')
    set(gca, 'FontSize', fontsize)
    title(param_update_names_wt{i_param}, 'FontSize', fontsize)
    grid on
end

%% ksdensity and prior
figure
set(gcf, 'position', [298 107 1259 822])
for i_param = 1:n_update_param_wt
    subplot(n_row, n_col, i_param)
    [f,xi] = ksdensity(log10(map_param_vals_wt(:, i_param)));
    plot(xi, f, '-', 'linewidth', 2)
    hold on
    title(param_update_names_wt{i_param})
    grid on
end

% show parameter prior distribution
load(sorted_wt{1,'filepath'}{1}, 'parameter_update');
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

%% show the steady state concentrations for 12 species and free+complex
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1', 'Mig1*', 'C83', 'C84', 'glu', 'gal'};
n_species = length(species_list);

figure
set(gcf, 'position', [298 107 1182 856])

for i_species = 1:n_species
    subplot(4, 4, i_species)
    plot([1:8], wt_conc_Gal_list.col(:,:,i_species), '.', 'markersize', markersize)
    hold on
    plot([1:12], wt_conc_Gal_list.row(:,:,i_species), '.', 'markersize', markersize)
    xlim([0 13])
%     set(gca, 'yscale', 'log')
    set(gca, 'FontSize', fontsize)
    title(species_list{i_species}, 'FontSize', fontsize)
    grid on
end

% plot free+complex steady state

subplot(4,4,13) % G3+G3*+C83
plot(1:8, (wt_conc_Gal_list.col(:,:,3)+wt_conc_Gal_list.col(:,:,6)+wt_conc_Gal_list.col(:,:,9)), '.', 'markersize', markersize)
hold on
plot(1:12, (wt_conc_Gal_list.row(:,:,3)+wt_conc_Gal_list.row(:,:,6)+wt_conc_Gal_list.row(:,:,9)), '.', 'markersize', markersize)
xlim([0 13])
set(gca, 'FontSize', fontsize)
title('G3+G3*+C83', 'FontSize', fontsize)
grid on

subplot(4,4,14) % G80+C83+C84
plot(1:8, (wt_conc_Gal_list.col(:,:,5)+wt_conc_Gal_list.col(:,:,9)+wt_conc_Gal_list.col(:,:,10)), '.', 'markersize', markersize)
hold on
plot(1:12, (wt_conc_Gal_list.row(:,:,5)+wt_conc_Gal_list.row(:,:,9)+wt_conc_Gal_list.row(:,:,10)), '.', 'markersize', markersize)
xlim([0 13])    
set(gca, 'FontSize', fontsize)
title('G80+C83+C84', 'FontSize', fontsize)
grid on

subplot(4,4,15) % G4+C84
plot(1:8, (wt_conc_Gal_list.col(:,:,4)+wt_conc_Gal_list.col(:,:,10)), '.', 'markersize', markersize)
hold on
plot(1:12, (wt_conc_Gal_list.row(:,:,4)+wt_conc_Gal_list.row(:,:,10)), '.', 'markersize', markersize)
xlim([0 13])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

subplot(4,4,16) % R+R*
plot(1:8, (wt_conc_Gal_list.col(:,:,7)+wt_conc_Gal_list.col(:,:,8)), '.', 'markersize', markersize)
hold on
plot(1:12, (wt_conc_Gal_list.row(:,:,7)+wt_conc_Gal_list.row(:,:,8)), '.', 'markersize', markersize)
xlim([0 13])
set(gca, 'FontSize', fontsize)
title('Mig1+Mig1*', 'FontSize', fontsize)
grid on

