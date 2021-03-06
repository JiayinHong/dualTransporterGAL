% this script is used to load and plot the results generated by function
% 'singleParamScan', namely obj versus multi-strain single param scan
% plotting, as well as visualization of 96 well trait using heatmap

%% load obj results from each file, then store them all in 'obj_all'
clear
obj_all = struct;
tmp = dir('../paramScanResults/');
for i=3:length(tmp)     % the first 2 are not files
    if ~ strcmp(tmp(i).name, '.DS_Store')
        load(fullfile('../paramScanResults/', tmp(i).name))
        str = fieldnames(obj);  % fieldnames is a cell
        str = str{1};           % convert cell to char
        no = str(4:end);        % 'str10' -> '10'
        %         no = str2num(no)-2;     % '10' -> 8
        no = str2num(no)-3;     % be careful to -2 / -3,
        % if there's '.DS_Store' in the traits folder, should -3
        fdn = sprintf('str%02d',no);    % fdn='str08'
        obj_all.(fdn) = obj.(str);  % assign the value to the same field of obj_all
    end
end

%% store all the strain names into a cell array
clear strNames
nStr = numel(fieldnames(obj_all));
strNames = cell(nStr,1);

for i=3:length(allTraits)   % the first 2 are not traits
    if ~ strcmp(allTraits(i).name, '.DS_Store')
        tmp = allTraits(i).name;
        %         strNames{i-2} = tmp(1:end-4);   % to remove '.mat'
        strNames{i-3} = tmp(1:end-4);   % to remove '.mat'
    end
end

%% load base parameters (best-fit) of wildtype
load('../wildtype_96well-025-170925_16:18.mat')
param_names = parameter_update.parameter_name;
base_param = param_map;

% strain sets for data from Kayla (multi strain comparison)
% strainSet1 = {'str13','str22','str08','str10','str16'};
% strainSet2 = {'str34','str23','str26','str27'};
% strainSet3 = {'str01','str02','str11','str12','str29'};
% strainSet4 = {'str32','str05','str20'};
% strainContainNaN = {'str04','str06','str19','str31','str33'};
% outpath = '../multiStrScan/';

% strain sets for simulation generated fake traits
for i=1:8
    MMconstAct{i} = sprintf('str%02d',i);
end
for i=9:16
    MMconstRep{i-8} = sprintf('str%02d',i);
end
for i=17:26
    basalSynthesis{i-16} = sprintf('str%02d',i);
end
for i=27:34
    induceSynthesis{i-26} = sprintf('str%02d',i);
end
for i=35:40
    Others{i-34} = sprintf('str%02d',i);
end
for i=41:48
    bindRate{i-40} = sprintf('str%02d',i);
end
for i=49:58
    hillCoeff{i-48} = sprintf('str%02d',i);
end

%% plot multi-strains single param scan on one fig
if 1
    outpath = '../fakeTraitsScan/';
    strainSet = hillCoeff;
    setName = 'hillCoeff';
    for i_param = 1:length(param_names)
        param_name = param_names{i_param};
        base_val = base_param.(param_name);
        multiStrParamScanHelper(param_name, strNames, base_val, obj_all, i_param, outpath, strainSet, setName)
    end
end

%% this section is used to visualize the trait in a heatmap
for i=3:length(allTraits)
    tmp = allTraits(i).name;
    filename = tmp(1:end-4);
    load(sprintf('../traitExtraction/%s', tmp))
    target_trait = trait;
    
    plot_heatmap(filename, target_trait)
end

%%  below is in-script functions

function multiStrParamScanHelper(param_name, strNames, base_val, obj_all, i_param, outpath, strainSet, setName)

obj_thresh = 10;    % to draw a reference line
cmap = lines(7);    % customize a colormap, the default colomap for lines is 'lines' which contains 7 unique colors
cmap(8:10,:) = winter(3);   % in our case, it's not sufficient, so I add three more colors from other colormaps
cmap(9,:) = autumn(1);  % for more details and options, see manual of 'colormap'

if regexp(param_name, 'n*')     % hill coefficients
    perturb_coefficient = linspace(1,4,7);
    nValue = length(perturb_coefficient);
    varied_value = perturb_coefficient;
    catnames = cell(1,nValue);
    for i=1:nValue
        catnames{i} = num2str(varied_value(i), '%.1f');
    end
    
else
    perturb_coefficient = logspace(log10(0.001), log10(1000), 13);
    nValue = length(perturb_coefficient);
    varied_value = perturb_coefficient .* base_val;
    catnames = cell(1,nValue);
    for i=1:nValue
        catnames{i} = num2str(varied_value(i), '%.4f');
    end
    
end

figure
set(gcf, 'position', [2049 418 1065 411])
whichStr=[];
j=1;
for fdn = strainSet                % select a couple of strains to draw on the same plot
    tmp = obj_all.(fdn{1});
    plot(1:nValue, tmp(i_param,1:nValue), 'Color', cmap(j,:), 'LineWidth', 1.5, 'Marker', 'o')
    whichStr(j) = str2num(fdn{1}(4:end));   % iStr stores the strains that have been drawn
    j=j+1;
    hold all
end

hline = refline([0 obj_thresh]);    % add a reference line
hline.LineWidth = 1.5;
hline.Color = 'm';
hline.LineStyle = '-.';

set(gca, 'xtick', 1:nValue, 'xticklabels', catnames)
if regexp(param_name, 'n*')     % hill coefficients
    % no rotation
else
    set(gca, 'XTickLabelRotation', 45)
end
grid on
legend(strNames(whichStr), 'location', 'NorthWest')
xlabel(sprintf('%s param values', param_name))
ylabel('sum obj - 96 well')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')

export_fig(fullfile(outpath, sprintf('%s vary in logspace - %s', param_name, setName)));
end

function plot_heatmap(filename, target_trait, param_set)
%   warning: be careful of the trait loading, there is only once load in
%   the current version, if changed later, be sure the old trait are not
%   overlaid by the new one

%% claim saving directory
saveDir = sprintf('../multiStrain96wellPlot/');  % the directory to store the figures
if ~isdir(saveDir)
    mkdir(saveDir)
end

%% load experimental data

expt_96well = target_trait;
galLabel = {'None','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2'};
gluLabel = {'None','-6','-5','-4','-3','-2','-1','0'};
colLabels = fliplr(gluLabel);
rowLabels = galLabel;
load_global
alldata = nan(8,12);

% to clean the data
ind1 = find(expt_96well.mask_induction == 0);     % all the rows whose mask_induction == 0
tmp = find(expt_96well(ind1,:).mask_basal == 0);
ind2 = ind1(tmp);                               % the rows whose mask_basal also equals to 0
ind1(tmp) = [];                                 % remove ind2 from ind1, so that ind1 only contains
% rows whose mask_induction == 0 while mask_basal ~= 0

% use NaN for those induced level in ind2
expt_96well{ind2, 'ind_level'} = NaN;
% use basal_level to represent induced level in ind1
expt_96well(ind1,:).ind_level = expt_96well(ind1,:).basal_level;

alldata(:) = logyfp_to_nm(expt_96well{:,'ind_level'});

%% heatmap for the expt trait
figure
set(gcf, 'position', [689 136 1036 811])
heatmap(rowLabels, colLabels, alldata, 'CellLabelFormat', '%.2f', 'FontSize', 12);
title(sprintf('%s expt G1 induced level', filename));
export_fig(fullfile(saveDir, get(gca, 'Title')))

if nargin==3    % there is a set of parameters that user wants to
    % pass to the function and calculate the difference
    
    %% fetch simulation results
    output = evalGalPathway( param_set, target_trait, '96well' );
    
    simG1_96well = nan(8,12);
    simG1_96well(:) = output.all_conc_Gal(:,1);
    
    if exist('ind1', 'var')
        % use basal level to compare with expt data when mask_basal = 1 &&
        % mask_induction == 0
        simG1_96well(ind1) = output.all_conc_Glu(ind1,1);
    end
    
    % get GAL1 obj
    GAL1_obj = output.G1obj;
    
    %% simple calculation
    logAllData = log10(alldata);
    logSimG1 = log10(simG1_96well);
    logdelta = logSimG1 - logAllData;   % the deviation in log scale
    fprintf('\nthe total difference: %.2f\n', nansum(abs(logdelta(:))))
    
    %% draw the other two heatmap
    
    % second, heatmap for simulation results
    figure
    set(gcf, 'position', [689 136 1036 811])
    heatmap(rowLabels, colLabels, simG1_96well, 'CellLabelFormat', '%.2f', 'FontSize', 12);
    title(sprintf('%s simulation G1 induced level', filename));
    export_fig(fullfile(saveDir, get(gca, 'Title')))
    
    % then, difference heatmap
    figure
    set(gcf, 'position', [689 136 1036 811])
    heatmap(rowLabels, colLabels, logdelta, 'ColorLimits', [-1 1] ...
        ,'CellLabelFormat', '%.2f', 'FontSize', 12);
    title(sprintf('The deviation heatmap, obj=%.2f', GAL1_obj));
    export_fig(fullfile(saveDir, get(gca, 'Title')))
    
end

end
