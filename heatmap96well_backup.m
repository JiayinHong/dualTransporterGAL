%% heatmap showing experimental wildtype 96-well induced G1 level

load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat')
G1_96well = trait;
galLabel = {'None','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2'};
gluLabel = {'None','-6','-5','-4','-3','-2','-1','0'};
colLabels = fliplr(gluLabel);
rowLabels = galLabel;
load_global
alldata = nan(8,12);
% to clean the data
ind1 = find(G1_96well.mask_induction == 0);     % all the rows whose mask_induction == 0
tmp = find(G1_96well(ind1,:).mask_basal == 0);
ind2 = ind1(tmp);                               % the rows whose mask_basal also equals to 0
ind1(tmp) = [];                                 % remove ind2 from ind1, so that ind1 only contains
% rows whose mask_induction == 0 while mask_basal ~= 0
% use mean to represent induced level in ind2
% G1_96well(ind2,:).ind_level = G1_96well(ind2,:).basal_level .* G1_96well(ind2,:).basal_frac ...
%     + G1_96well(ind2,:).ind_level .* G1_96well(ind2,:).ind_frac;

% use NaN for those induced level in ind2
G1_96well{ind2, 'ind_level'} = NaN;
% use basal_level to represent induced level in ind1
G1_96well(ind1,:).ind_level = G1_96well(ind1,:).basal_level;

alldata(:) = logyfp_to_nm(G1_96well{:,'ind_level'});


figure
set(gcf, 'position', [689 136 1036 811])
% R2016a version
% heatmap(alldata, rowLabels, colLabels, '%.2f', 'Colorbar', true ...
%     , 'ShowAllTicks', true, 'TextColor', 'r', 'FontSize', 14  ...
%     , 'GridLines', ':', 'ColorLevels', 128, 'TickFontSize', 15);
% title('Wildtype G1 induced level', 'FontSize', 15)
% xlabel('galactose titration')
% ylabel('glucose titration')

% R2017a version
h = heatmap(rowLabels, colLabels, alldata, 'CellLabelFormat', '%.2f', 'FontSize', 15);
title('Wildtype G1 induced level')

%% Load mcmc results
load('../metaData/trait_extraction/GAL3pr_all_data.mat')
G3level = trait.basal_level;
load('../metaData/trait_extraction/GAL4pr_all_data.mat')
G4level = trait.basal_level;

mcmc_data_folder = '../results/biTrans_addHXT_96well/';
jobtags = {'wildtype_96well'};
biTrans_addHXT_96well = load_mcmc_result(mcmc_data_folder, jobtags);

%% heatmap showing simulated wildtype 96-well induced G1 level
i_example = 1;
param = biTrans_addHXT_96well{i_example, 'param_map'};
output = evalGalPathway_GAL34_changedR(param, G1_96well, 0, G3level, G4level, '96well');

simG1_96well = nan(8,12);
simG1_96well(:) = output.all_conc_Gal(:,1);
% use basal level to compare with expt data when mask_basal = 1 &&
% mask_induction == 0
simG1_96well(ind1) = output.all_conc_Glu(ind1,1);
% fetch original simulation result to draw the titration plot
simG1_ind = output.all_conc_Gal(:,1);
simG1_basal = output.all_conc_Glu(:,1);

%%
figure
set(gcf, 'position', [689 136 1036 811])
% R2016a version

% heatmap(simG1_96well, rowLabels, colLabels, '%.2f' ...
%     , 'ShowAllTicks', true, 'TextColor', 'r', 'FontSize', 14, 'Colorbar', true ...
%     , 'GridLines', ':', 'ColorLevels', 128, 'TickFontSize', 15);
% title(sprintf('The no.%d example', i_example), 'FontSize', 15)
% xlabel('galactose titration')
% ylabel('glucose titration')

% R2017a version
heatmap(rowLabels, colLabels, simG1_96well, 'CellLabelFormat', '%.2f', 'FontSize', 15);
title(sprintf('The no.%d example', i_example))


%% difference map
% experimental data 8*12 - alldata
% simulation results 8*12 - simG1_96well
logAllData = log(alldata);
logSimG1 = log(simG1_96well);
logdelta = logSimG1 - logAllData;   % the deviation in log scale
lindelta = simG1_96well - alldata;  % the deviation in linear scale

figure
set(gcf, 'position', [689 136 1036 811])

% R2016a version
% heatmap(logdelta, rowLabels, colLabels, '%.2f' ...
%     , 'ShowAllTicks', true, 'TextColor', 'r', 'FontSize', 14, 'Colorbar', true ...
%     , 'GridLines', ':', 'ColorLevels', 128, 'TickFontSize', 15);
% title(sprintf('The deviation heatmap of no.%d example', i_example), 'FontSize', 15)
% xlabel('galactose titration')
% ylabel('glucose titration')

% R2017a version
heatmap(rowLabels, colLabels, logdelta, 'CellLabelFormat', '%.2f', 'FontSize', 15);
title(sprintf('The deviation heatmap of no.%d example', i_example))


%% compress the data in either direction
% reload wildtype 96-well G1 level data
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat')
markersize = 6;
linewid = 1.5;

% split the 96-well plate into 8 rows, each one is galactose titration
galTitrate = 1:8:89;    % the first row of gal titration

figure
set(gcf, 'position', [680 106 1052 872])

for i = 1:8     % from the first to the last row
    subplot(8,1,i)
    for j = 1:12
        index = i - 1 + galTitrate;   % get the subscript in the trait table
        sub = index(j);
        % plot experimental basal
        if trait{sub, 'mask_basal'}
            plot(j,logyfp_to_nm(trait{sub, 'basal_level'}), 'ok', 'markersize', markersize)
        else
            plot(j,logyfp_to_nm(trait{sub, 'basal_level'}), '+k', 'markersize', markersize)
        end
        hold all
        % plot experimental induced
        if trait{sub, 'mask_induction'}
            plot(j,logyfp_to_nm(trait{sub, 'ind_level'}), 'or', 'markersize', markersize)
        else
            plot(j,logyfp_to_nm(trait{sub, 'ind_level'}), '+r', 'markersize', markersize)
        end
    end
    
    plot(simG1_basal(index), 'k-', 'linewidth', linewid)
    plot(simG1_96well(index), 'r-', 'linewidth', linewid)
    
    set(gca, 'yscale', 'log', 'FontSize', 12)
    if i == 8       % the last row
        set(gca, 'XTick', 0:13)
        set(gca, 'XTickLabel', {'', galLabel{:}, ''})
    else
        set(gca, 'XTickLabel', [])
    end
    xlim([0, 13])
    ylabel(colLabels{i}, 'FontWeight', 'bold')
end

h = suplabel('galactose titration subplot', 't');
h.FontSize = 15;
h = suplabel('glucose gradient', 'y');
h.FontSize = 15;
h = suplabel('galactose gradient', 'x');
h.FontSize = 15;


%% split the 96-well plate into 12 rows, each one is glucose titration
gluTitrate = 1:8;    % the first col of glu titration

figure
set(gcf, 'position', [680 106 1052 872])

for i = 1:12     % from the first to the last row
    subplot(12,1,i)
    for j = 1:8
        index = gluTitrate + 8*(i-1);   % get the subscript in the trait table
        sub = index(j);
        % plot experimental basal
        if trait{sub, 'mask_basal'}
            plot(j,logyfp_to_nm(trait{sub, 'basal_level'}), 'ok', 'markersize', markersize)
        else
            plot(j,logyfp_to_nm(trait{sub, 'basal_level'}), '+k', 'markersize', markersize)
        end
        hold all
        % plot experimental induced
        if trait{sub, 'mask_induction'}
            plot(j,logyfp_to_nm(trait{sub, 'ind_level'}), 'or', 'markersize', markersize)
        else
            plot(j,logyfp_to_nm(trait{sub, 'ind_level'}), '+r', 'markersize', markersize)
        end
    end
    
    plot(simG1_basal(index), 'k-', 'linewidth', linewid)
    plot(simG1_96well(index), 'r-', 'linewidth', linewid)
    
    set(gca, 'yscale', 'log', 'FontSize', 12)
    if i == 12      % the last row
        set(gca, 'XTickLabel', {'',colLabels{:},''})
    else
        set(gca, 'XTickLabel', [])
    end
    xlim([0, 9])
    ylabel(rowLabels{i}, 'FontWeight', 'bold')
end

h = suplabel('glucose titration subplot', 't');
h.FontSize = 15;
h = suplabel('galactose gradient', 'y');
h.FontSize = 15;
h = suplabel('glucose gradient', 'x');
h.FontSize = 15;


%% split the 96-well plate into 12 columns, each one is glucose titration
gluTitrate = 1:8;    % the first col of glu titration

figure
set(gcf, 'position', [446 106 1286 872])

for i = 1:12     % from the first to the last row
    subplot(1,12,i)
    for j = 1:8
        index = gluTitrate + 8*(i-1);   % get the subscript in the trait table
        sub = index(j);
        % plot experimental basal
        if trait{sub, 'mask_basal'}
            plot(logyfp_to_nm(trait{sub, 'basal_level'}),9-j, 'ok', 'markersize', markersize)
        else
            plot(logyfp_to_nm(trait{sub, 'basal_level'}),9-j, '+k', 'markersize', markersize)
        end
        hold all
        % plot experimental induced
        if trait{sub, 'mask_induction'}
            plot(logyfp_to_nm(trait{sub, 'ind_level'}),9-j, 'or', 'markersize', markersize)
        else
            plot(logyfp_to_nm(trait{sub, 'ind_level'}),9-j, '+r', 'markersize', markersize)
        end
    end
    
    plot(simG1_basal(index),8:-1:1, 'k-', 'linewidth', linewid)
    plot(simG1_96well(index),8:-1:1, 'r-', 'linewidth', linewid)
    
    set(gca, 'xscale', 'log', 'FontSize', 12)
    if i == 1   % the first col
        set(gca, 'YTickLabel', {'','None','-6','-5','-4','-3','-2','-1','0',''})
    else
        set(gca, 'YTickLabel', [])
    end
    ylim([0, 9])
    xlabel(rowLabels{i}, 'FontWeight', 'bold')
end

h = suplabel('glucose titration subplot', 't');
h.FontSize = 15;
h = suplabel('galactose gradient', 'x');
h.FontSize = 15;
h = suplabel('glucose gradient', 'y');
h.FontSize = 15;




