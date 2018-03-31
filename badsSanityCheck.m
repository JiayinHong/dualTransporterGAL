%% load bads results
if 1
    badsResult = read_bads_result('../results/badsOptim/RenanData/');
    badsResult = sortrows(badsResult,'sum_obj','ascend');
    
    % load configurations, i.e. parameter_update, param_init, etc
    folder_name = '../metaData/noCompeteBinding/';
    jobtag = 'BC187_Kayla';     % choose from 'BC187_Kayla', 'BC187_Renan' & 'YJM978_Renan'
    filepath = fullfile(folder_name, [jobtag, sprintf('_%03d',str2num('1')), '.mat']);  % config files are identical for the same batch
    load(filepath)
end

%% plot 96 well heatmap for sanity check
i_example = 11;
parameter_name = parameter_update.parameter_name;
% call function 'update_param' to get the full struct of parameter values
best_param = update_param( param_init, parameter_name, badsResult{i_example,'best_params'}{1});
% plot the heatmap
% heatmap96well(best_param,'bc187R','R2016a')     % DON'T forget to change the options here!

% save the best params and param names into a .mat file
parameter_update = readtable('Nov15th_param_config_set11.csv');
param_names = parameter_update.parameter_name;
tok = regexp(badsResult{i_example,'filename'}{1}, '([\w]+-\d{3})-*', 'tokens');
tok = tok{1};
save(['../bestParams/',tok{1},'.mat'],'best_param','param_names')

%% Steady State plot
% load GAL3, GAL4 level
load('../metaData/trait_extraction/GAL3pr_all_data.mat')
G3level = trait.basal_level;
load('../metaData/trait_extraction/GAL4pr_all_data.mat')
G4level = trait.basal_level;

% load parameter, trait, and fit type
jobtag = badsResult{i_example, 'jobtag'};
jobtag = jobtag{1};
param = best_param;
filepath = fullfile(folder_name, [jobtag, sprintf('_%03d',str2num('1')), '.mat']);
load(filepath)      % load trait
GAL1_trait = trait;
fit_type = 'Rcross';

% plot for 11 species excluding R*

clear lo_state hi_state Mig1star
figure
set(gcf, 'position', [298 107 1182 856])
markersize = 10;
fontsize = 12;

subplot(4,5,1)
% show how GAL1 level fits, also get steady state concentrations for 11
% species
[lo_state(:,:), hi_state(:,:)] = param_fitting_plot(param, GAL1_trait, fit_type);
subplot(4,5,6)
% show how GAL3 level fits
param_fitting_plot_GAL234(param, GAL1_trait, 0, G3level, G4level, fit_type, 'G3');
subplot(4,5,11)
% show how GAL4 level fits
param_fitting_plot_GAL234(param, GAL1_trait, 0, G3level, G4level, fit_type, 'G4');

% show the steady state concentrations for 11 species and free+complex
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1tot', 'C83', 'C84', 'glu', 'gal'};
n_species = length(species_list);
switch fit_type
    case 'one_row'
        n_condition = 12;
    case 'one_column'
        n_condition = 8;
    case 'one_cross'
        n_condition = 19;
    case 'Rcross'
        n_condition = 18;
end

subplot(4,5,2)
plot(1:n_condition, hi_state(:,1), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
title('G1', 'FontSize', fontsize)
grid on

ax9 = subplot(4,5,9);   % return the axis object of subplot(4,5,9),
                        % the one of species 'Mig1tot'

for i_species = 2:n_species
    if i_species > 8
        subplot(4, 5, i_species+3)
    elseif i_species > 4
        subplot(4, 5, i_species+2)
    else
        subplot(4, 5, i_species+1)
    end
    
    plot(1:n_condition, hi_state(:,i_species), '.', 'markersize', markersize)
    
    xlim([0 n_condition+1])
    %     set(gca, 'yscale', 'log')
    set(gca, 'FontSize', fontsize)
    title(species_list{i_species}, 'FontSize', fontsize)
    grid on
end

% ytickformat(ax9,'%.0f')     % specify the tick label format of 'Mig1tot'

% in matlab R2016a, the code above is not functional, use code below
shortTickLabel = cellfun(@(s) sprintf('%.6s',s), ax9.YTickLabel, 'UniformOutput',false);
ax9.YTickLabel = shortTickLabel;

% plot Mig1* level
subplot(4,5,15)
Mig1tot = hi_state(:,7);
glu_in = hi_state(:,10);
if range(Mig1tot) == 0 || range(Mig1tot) < 1e-7  % steady state level = aR/gamma, no matter what condition
    Mig1star = glu_in .^ param.nRs ./ (param.KRs^param.nRs + glu_in .^ param.nRs) .* Mig1tot(1);
else
    error('There''s something wrong!')
end
plot(1:n_condition, Mig1star, '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('Mig1*', 'FontSize', fontsize)
grid on

% plot free+complex steady state

subplot(4,5,17) % G3+G3*+C83

plot(1:n_condition, (hi_state(:,3) + hi_state(:,6) + hi_state(:,8)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G3+G3*+C83', 'FontSize', fontsize)
grid on

subplot(4,5,18) % G80+C83+C84

plot(1:n_condition, (hi_state(:,5) + hi_state(:,8) + hi_state(:,9)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G80+C83+C84', 'FontSize', fontsize)
grid on

subplot(4,5,19) % G4+C84

plot(1:n_condition, (hi_state(:,4) + hi_state(:,9)), '.', 'markersize', markersize)
% ytickformat('%.2f')     % specify the tick label format of 'G4+C84'
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

jobtag = changeunderscore(jobtag);
[ax,h] = suplabel(sprintf('The no.%s example, %s', num2str(i_example), jobtag), 't');
h.FontSize = 13;


