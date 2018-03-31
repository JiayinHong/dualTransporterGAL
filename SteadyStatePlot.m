%% load GAL3, GAL4 level
load('../metaData/trait_extraction/GAL3pr_all_data.mat')
G3level = trait.basal_level;
load('../metaData/trait_extraction/GAL4pr_all_data.mat')
G4level = trait.basal_level;

%%
% mcmc_data_folder = '../results/biTrans_addHXT/';
% jobtags = {'medium-wildtype_1c', 'medium-wildtype_1r', 'medium-wildtype_1r1c'};
% mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);


%% load parameter, trait, and fit type

i_example = 1;
jobtag = mcmc_result{i_example, 'jobtag'};
param = mcmc_result{i_example, 'param_map'};
% trait = mcmc_result{i_example, 'trait'}{1};
load(mcmc_result{i_example,'filepath'}{1},'GAL1_trait')
fit_type = 'one_cross';

%% plot for 12 species including HXT
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
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1tot', 'C83', 'C84', 'glu', 'gal', 'HXT'};
n_species = length(species_list);
switch fit_type
    case 'one_row'
        n_condition = 12;
    case 'one_column'
        n_condition = 8;
    case 'one_cross'
        n_condition = 19;
end

subplot(4,5,2)
plot(1:n_condition, hi_state(:,1), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
title('G1', 'FontSize', fontsize)
grid on

ax15 = subplot(4,5,15);

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

ax15.YScale = 'log';

% plot Mig1* level
subplot(4,5,20)
Mig1tot = hi_state(:,7);
glu_in = hi_state(:,10);
if range(Mig1tot) == 0  % steady state level = aR/gamma, no matter what condition
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
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

jobtag = changeunderscore(jobtag{1});
[ax,h] = suplabel(sprintf('The no.%s example, %s', num2str(i_example), jobtag), 't');
h.FontSize = 13;

saveDir = '../steadyStatePlot/';
if ~isdir(saveDir)
    mkdir(saveDir)
end
% export_fig(fullfile(saveDir, h.String))
export_fig(fullfile(saveDir, 'wildtype96well'))

%% load parameter, trait, and fit type

i_example = 2;
jobtag = mcmc_old_result{i_example, 'jobtag'};
param = mcmc_old_result{i_example, 'param_map'};
trait = mcmc_old_result{i_example, 'trait'}{1};
fit_type = 'one_column';

%% plot for 11 species excluding R*
% remove Mig1* feedback in the system, glucose direct inhibit GAL pathway
clear lo_state hi_state Mig1star
figure
set(gcf, 'position', [298 107 1182 856])
markersize = 10;
fontsize = 12;

subplot(4,5,1)
% show how GAL1 level fits, also get steady state concentrations for 11
% species
[lo_state(:,:), hi_state(:,:)] = param_fitting_plot(param, trait, fit_type);
subplot(4,5,6)
% show how GAL3 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G3');
subplot(4,5,11)
% show how GAL4 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G4');

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
end

subplot(4,5,2)
plot(1:n_condition, hi_state(:,1), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
title('G1', 'FontSize', fontsize)
grid on

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
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

jobtag = changeunderscore(jobtag{1});
[ax,h] = suplabel(sprintf('The no.%s example, %s', num2str(i_example), jobtag), 't');
h.FontSize = 13;

%% plot for 11 species excluding R*
clear lo_state hi_state Mig1star
figure
set(gcf, 'position', [298 107 1182 856])
markersize = 10;
fontsize = 12;

subplot(4,5,1)
% show how GAL1 level fits, also get steady state concentrations for 11
% species
[lo_state(:,:), hi_state(:,:)] = param_fitting_plot(param, trait, fit_type);
subplot(4,5,6)
% show how GAL3 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G3');
subplot(4,5,11)
% show how GAL4 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G4');

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
end

subplot(4,5,2)
plot(1:n_condition, hi_state(:,1), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
title('G1', 'FontSize', fontsize)
grid on

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

% plot Mig1* level
subplot(4,5,15)
Mig1tot = hi_state(:,7);
glu_in = hi_state(:,10);
if range(Mig1tot) == 0  % steady state level = aR/gamma, no matter what condition
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
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

jobtag = changeunderscore(jobtag{1});
[ax,h] = suplabel(sprintf('The no.%s example, %s', num2str(i_example), jobtag), 't');
h.FontSize = 13;

%% plot for 12 species including R*
clear lo_state hi_state
figure
set(gcf, 'position', [298 107 1182 856])
markersize = 10;
fontsize = 12;

subplot(4,5,1)
% show how GAL1 level fits, also get steady state concentrations for 12
% species
[lo_state(:,:), hi_state(:,:)] = param_fitting_plot(param, trait, fit_type);
subplot(4,5,6)
% show how GAL3 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G3');
subplot(4,5,11)
% show how GAL4 level fits
param_fitting_plot_GAL234(param, trait, 0, G3level, G4level, fit_type, 'G4');

% show the steady state concentrations for 12 species and free+complex
species_list = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'Mig1', 'Mig1*', 'C83', 'C84', 'glu', 'gal'};
n_species = length(species_list);
switch fit_type
    case 'one_row'
        n_condition = 12;
    case 'one_column'
        n_condition = 8;
    case 'one_cross'
        n_condition = 19;
end


for i_species = 1:n_species
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

% plot free+complex steady state

subplot(4,5,17) % G3+G3*+C83

plot(1:n_condition, (hi_state(:,3) + hi_state(:,6) + hi_state(:,9)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G3+G3*+C83', 'FontSize', fontsize)
grid on

subplot(4,5,18) % G80+C83+C84

plot(1:n_condition, (hi_state(:,5) + hi_state(:,9) + hi_state(:,10)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G80+C83+C84', 'FontSize', fontsize)
grid on

subplot(4,5,19) % G4+C84

plot(1:n_condition, (hi_state(:,4) + hi_state(:,10)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('G4+C84', 'FontSize', fontsize)
grid on

subplot(4,5,20) % R+R*

plot(1:n_condition, (hi_state(:,7) + hi_state(:,8)), '.', 'markersize', markersize)
xlim([0 n_condition+1])
set(gca, 'FontSize', fontsize)
title('Mig1+Mig1*', 'FontSize', fontsize)
grid on

