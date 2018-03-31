%% load mcmc results
mcmc_data_folder = '../results/biTrans_addHXT_96well/';
% mcmc_data_folder = '../results/fixOverlap/';
jobtags = {'wildtype_96well', 'gal80d_96well', 'mig1d_96well'};
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);
mcmc_result = sortrows(mcmc_result, 'map_data_over_param', 'descend');

%% decide a threshold to filter data
prob_thresh = -50;
mcmc_result_filter = mcmc_result(mcmc_result.map_data_over_param > prob_thresh,:);

mcmc_result = sortrows(mcmc_result, 'map_data_over_param', 'descend');
mcmc_wt = mcmc_result(ismember(mcmc_result.jobtag, 'wildtype_96well'),:);
mcmc_mig1d = mcmc_result(ismember(mcmc_result.jobtag, 'mig1d_96well'),:);
mcmc_gal80d = mcmc_result(ismember(mcmc_result.jobtag, 'gal80d_96well'),:);
mcmc_result_filter = [mcmc_wt(1:16,:);mcmc_gal80d(1:16,:);mcmc_mig1d(1:16,:)];

%% fetch parameter values - only MAP
% we are now analyzing all the data that are above a probability threshold,
% not just MAP, so do not run this hunk of code for now.
if 0
    param_val_all = [];
    param_id_all = [];
    
    for i_job = 1:length(jobtags)
        
        mcmc_result_tmp = mcmc_result_filter(ismember(mcmc_result_filter.jobtag, jobtags{i_job}),:);
        load(mcmc_result_tmp{1, 'filepath'}{1}, 'parameter_update');
        
        parameter_update_names = parameter_update.parameter_name;
        n_sample = height(mcmc_result_tmp);
        n_param_update = height(parameter_update);
        param_map_val = nan(n_sample, n_param_update);
        
        for i_sample = 1:n_sample
            param_map = mcmc_result_tmp{i_sample, 'param_map'};
            for i_param_update = 1:n_param_update
                param_map_val(i_sample, i_param_update) = log10(param_map.(parameter_update_names{i_param_update}));
            end
        end
        
        % use id to distinguish between wt and mutants clusters
        param_map_id = repmat(i_job,[n_sample,1]);
        param_val_all = [param_val_all; param_map_val];
        param_id_all = [param_id_all; param_map_id];
    end
end

%% fetch parameter values - posterior distribution beyond a threshold
param_val_all = [];
param_id_all = [];
load(mcmc_result_filter{1, 'filepath'}{1}, 'parameter_update');
parameter_update(2,:)=[];

parameter_update_names = parameter_update.parameter_name;
n_param_update = height(parameter_update);

for i_job = 1:length(jobtags)
    
    mcmc_result_tmp = mcmc_result_filter(ismember(mcmc_result_filter.jobtag, jobtags{i_job}),:);
    n_chain = height(mcmc_result_tmp);      % the total Markov chains
    
    for i_chain = 1:n_chain
        prob_list = mcmc_result_tmp{i_chain, 'prob_data_over_parameter_list'}{1};
        param_list = mcmc_result_tmp{i_chain, 'param_list'};
        filter_id = find(prob_list > prob_thresh);
        [~,I] = max(prob_list);
        II = find(I==filter_id);    % the index of MAP
        param_chain = nan(length(filter_id), n_param_update);   % only contain the current strain values
        
        for i_param_update = 1:n_param_update
            tmp = param_list.(parameter_update_names{i_param_update});
            param_chain(:, i_param_update) = log10(tmp(filter_id));
        end
        param_val_all =  [param_val_all; param_chain];  % store all strains parameter values here
        param_chain_id = repmat(i_job, [length(filter_id),1]);
        param_chain_id(II) = 666;   % the ID of MAP
        param_id_all = [param_id_all; param_chain_id];  % store all ids for each set of parameters
    end
    
    
end


%% fetch parameter values - posterior distribution beyond a threshold
% different filter from last session
param_val_all = [];
param_id_all = [];
load(mcmc_result_filter{1, 'filepath'}{1}, 'parameter_update');
parameter_update_names = parameter_update.parameter_name;
n_param_update = height(parameter_update);

n_dataPoint = 500;

for i_job = 1:length(jobtags)
    
    mcmc_result_tmp = mcmc_result_filter(ismember(mcmc_result_filter.jobtag, jobtags{i_job}),:);
    n_chain = height(mcmc_result_tmp);      % the total Markov chains
    
    for i_chain = 1:n_chain
        prob_list = mcmc_result_tmp{i_chain, 'prob_data_over_parameter_list'}{1};
        param_list = mcmc_result_tmp{i_chain, 'param_list'};
        prob_list(isnan(prob_list)) = [];   % remove NaN in prob_list
        [~,filter_id] = sort(prob_list,'descend');
        
        param_chain = nan(n_dataPoint, n_param_update);   % only contain the current strain values
        
        for i_param_update = 1:n_param_update
            tmp = param_list.(parameter_update_names{i_param_update});
            param_chain(:, i_param_update) = log10(tmp(filter_id(1:n_dataPoint)));
        end
        param_val_all =  [param_val_all; param_chain];  % store all strains parameter values here
        param_chain_id = repmat(i_job, [n_dataPoint,1]);
        param_chain_id(1) = 666;   % the ID of MAP
        param_id_all = [param_id_all; param_chain_id];  % store all ids for each set of parameters
    end
    
    
end


%% remove the repetitions in param_val_all
[param_val_unique,ia,~] = unique(param_val_all, 'rows');   % param_val_unique contains all the unique values in the array
% now re-map the strain(cluster) ID to the new unique array
param_id_unique = param_id_all(ia);
% clear the redundant dataset
clear param_val_all
clear param_id_all

%% PCA and plot
[coeff, score, latent] = pca(zscore(param_val_unique));

fontsize = 12;
markersize = 15;

figure
plot([latent(1:10)/sum(latent) cumsum(latent(1:10))/sum(latent)]*100, 'LineWidth', 1.5);
set(gca,'FontSize', fontsize, 'FontWeight', 'bold')
xlabel('# of principal component');
ylabel('% of variance of dataset explained');
grid on
h=legend('Individual', 'Cumulative');
h.FontSize = fontsize;
h.FontWeight = 'bold';
h.Location = 'NorthWest';

%% 3d PCA
cluster_names = {'wt', 'gal80d', 'mig1d'};
figure
set(gcf,'position',[945 424 762 502])
for i_job = 1:length(jobtags)
    scatter3(score(param_id_unique==i_job,1),score(param_id_unique==i_job,2),score(param_id_unique==i_job,3),'.');
    hold on
%     text(mean(score(param_id_all == i_job,1)), mean(score(param_id_all==i_job,2)), mean(score(param_id_all == i_job,3)), cluster_names{i_job}, 'color', 'k', 'fontsize', 15, 'fontweight', 'bold');
    
end
set(gca,'FontSize', fontsize, 'FontWeight', 'bold')

% get the default RGB vector for each phenotype
h = get(gca);
wtC = h.ColorOrder(1,:);        % the color for wt
gal80dC = h.ColorOrder(2,:);    % the color for gal80d
mig1dC = h.ColorOrder(3,:);     % the color for mig1d

xlabel('1st PC')
ylabel('2nd PC')
zlabel('3rd PC')
legend(cluster_names, 'FontSize', fontsize, 'FontWeight', 'bold')
grid on
box on

%% 2d PCA highlight MAP
figure
set(gcf,'position',[484 152 1223 774])

subplot(2,2,1)
plot(score(:,1),score(:,2),'.','color',[0.7 0.7 0.7]);
hold on
plot(score(param_id_unique==666,1),score(param_id_unique==666,2),'d','MarkerFaceColor','red','markersize',3);
title('PC1 vs PC2', 'fontsize', 12)
legend({'non-MAP','MAP'}, 'FontSize', fontsize, 'FontWeight', 'bold', 'location', 'best')
grid on
set(gca,'fontsize',12)

subplot(2,2,2)
plot(score(:,1),score(:,3),'.','color',[0.7 0.7 0.7]);
hold on
plot(score(param_id_unique==666,1),score(param_id_unique==666,3),'d','MarkerFaceColor','red','markersize',3);
title('PC1 vs PC3', 'fontsize', 12)
legend({'non-MAP','MAP'}, 'FontSize', fontsize, 'FontWeight', 'bold', 'location', 'best')
grid on
set(gca,'fontsize',12)

subplot(2,2,3)
plot(score(:,2),score(:,3),'.','color',[0.7 0.7 0.7]);
hold on
plot(score(param_id_unique==666,2),score(param_id_unique==666,3),'d','MarkerFaceColor','red','markersize',3);
title('PC2 vs PC3', 'fontsize', 12)
legend({'non-MAP','MAP'}, 'FontSize', fontsize, 'FontWeight', 'bold', 'location', 'best')
grid on
set(gca,'fontsize',12)

% set(gca,'FontSize', fontsize, 'FontWeight', 'bold')

%% 2d PCA show density
pairscatter(score(:,1:3),[],[],'showDensity',true)
set(gcf, 'position', [360 130 780 568])

%% 2d PCA show cluster
figure
set(gcf,'position',[360 130 780 568])
for i_job = 1:length(jobtags)
    subplot(2,2,1)
    scatter(score(param_id_unique==i_job,1),score(param_id_unique==i_job,2),'.')
    title('PC1 vs PC2', 'fontsize', 12)
    set(gca,'fontsize',12)
    hold on
    
    subplot(2,2,2)
    scatter(score(param_id_unique==i_job,1),score(param_id_unique==i_job,3),'.')
    title('PC1 vs PC3', 'fontsize', 12)
    set(gca,'fontsize',12)
    hold on
    
    subplot(2,2,3)
    scatter(score(param_id_unique==i_job,2),score(param_id_unique==i_job,3),'.')
    title('PC2 vs PC3', 'fontsize', 12)
    set(gca,'fontsize',12)
    hold on
end
    
%% contribution of variables to Dim 1-nPC
nPC = 2;
nVar = length(parameter_update_names);

latentWeight = latent(:)/sum(latent);   % the weight of each PC
cutoff = sum(latentWeight(1:nPC))/nVar; % the average contribution as the cutoff
contrib = zeros(nVar,1);
for iPC = 1:nPC
    contribPerDim = coeff(:, iPC) .^2 * latentWeight(iPC);
    contrib = contrib+contribPerDim;    % weighted contribution to PC1-PCn of each original variable
end

figure
set(gcf, 'position', [360 120 777 578])
bar(contrib)
hold on
addzeroline2('ypos', cutoff, 'plotoption', {'color','r'})
set(gca, 'XTick', 1:nVar)
set(gca, 'XTickLabel', parameter_update_names)
set(gca, 'XTickLabelRotation', 45)
set(gca, 'FontSize', fontsize)
title(sprintf('Contribution of variables to Dim 1-%s', num2str(nPC)))

importantVar = find(contrib > cutoff);  % all the variables that possess a contribution above average
[~,topVar] = sort(contrib, 'descend');  % the index of each variable in a descending sequence of contribution

%% View the original variables (parameters that varied during MCMC)
% in the space of the first two or three principal components
figure
set(gcf, 'position', [360 90 849 608])
% the data points (observations) can also be shown in the same plot
% due to our vast dataset, omit the data points here
% for more information, check the documentation of biplot
h = biplot(coeff(:,1:2),  'varlabels', parameter_update_names); % the first two PC plane
set(gca, 'LineWidth', 1.5, 'FontSize', fontsize)

% the handle of biplot consists of 3*nVar+1 rows
% the first nVar rows are the lines, the second nVar rows are the points
% (the end of each line), the third nVar rows are texts, and the last one
% is the axis line.

% set all the lines wider than default
for i = 1:length(h)
    h(i,:).LineWidth = 1.5;
    % if there is a property called 'FontSize'
    if isprop(h(i,:), 'FontSize')
        h(i,:).FontSize = 12;
        h(i,:).FontWeight = 'bold';
    end
end

% if there are too many variables contributed above average
% only show top 10 varibales
if length(importantVar) > 10
    top10 = topVar(1:10);
    for i = [top10', nVar+top10']   % change the color of both lines and end points of lines
        h(i,:).Color = [1 0 0];
    end
    
else
    for i = [importantVar', nVar+importantVar'] % change the color of both lines and end points of lines
        h(i,:).Color = [1 0 0];
    end  
end

% change the width of axis line
h(end,:).LineWidth = 1;



%% violin plot revealing parameter value difference underlying wt and mutants
% only valid for plotting wt/mig1d/gal80d data

% first, get the subset parameter values from 'param_val_unique'
wtParamVal = param_val_unique(param_id_unique==1,:);
gal80dParamVal = param_val_unique(param_id_unique==2,:);
mig1dParamVal = param_val_unique(param_id_unique==3,:);

% then, initialize the structs
wtViolin = struct();
gal80dViolin = struct();
mig1dViolin = struct();

% store the value into the struct based on the fieldnames
for i_param_update = 1:n_param_update
    field_name = parameter_update_names{i_param_update};
    wtViolin.(field_name) = wtParamVal(:,i_param_update);
    gal80dViolin.(field_name) = gal80dParamVal(:,i_param_update);
    mig1dViolin.(field_name) = mig1dParamVal(:,i_param_update);
end

%% to free some memory
% clear wtParamVal
% clear gal80dParamVal
% clear mig1dParamVal

% now, make the plot
figure
set(gcf, 'position', [1946 145 1653 768])
hold all

%% Violin plot, usually takes ~20s
% It would be extremely slow if turning on the option - ShowData
% 'ViolinAlpha' is used to tune the transparency of the violins
% I hacked the function 'violinplot' to add an option 'Strain' in order to
% specify the plotting position for certain strain
% the color for all the violins should match the PCA plot

tic
violinplot(wtViolin,[], 'ViolinColor',wtC, 'Strain','wt', 'ShowData',false, 'ViolinAlpha',1);    % [ ] - no categories

violinplot(gal80dViolin,[], 'ViolinColor',gal80dC, 'Strain','gal80d', 'ShowData',false, 'ViolinAlpha',1);

[~,catnames] = violinplot(mig1dViolin,[], 'ViolinColor',mig1dC, 'Strain','mig1d', 'ShowData',false, 'ViolinAlpha',1);
toc

%% add a piece of code to extract the mean and median of parameters we're concerned with
% 2018.1.28 updated by JH
% pplist = {'a1','a4','a0HXT','ag4','d','kf3','KG1','KG80','KR3','KR4','KHXT','kG2','n1','n2','n3','n80','nRs'}';

% 2018.3.3 updated by JH
pplist = {'a2','a3','a80','aR','ag2','ag3','ag80','aHXT','KG2','KG3','KRs','rHXT','n1','nR1','nHXT'}';
pplist = {'ag1','kf83','kf84','KR1','rcat','rG2','KGglu','KHXTglu'}';

npp = length(pplist);

wt_mean = nan(npp,1);
wt_median = nan(npp,1);
gal80d_mean = nan(npp,1);
gal80d_median = nan(npp,1);
mig1d_mean = nan(npp,1);
mig1d_median = nan(npp,1);

for i=1:npp
    wt_mean(i) = 10^(mean(wtViolin.(pplist{i})));
    wt_median(i) = 10^(median(wtViolin.(pplist{i})));
    gal80d_mean(i) = 10^(mean(gal80dViolin.(pplist{i})));
    gal80d_median(i) = 10^(median(gal80dViolin.(pplist{i})));
    mig1d_mean(i) = 10^(mean(mig1dViolin.(pplist{i})));
    mig1d_median(i) = 10^(median(mig1dViolin.(pplist{i})));
end

param_summary = table(pplist,wt_mean,gal80d_mean,mig1d_mean,wt_median,gal80d_median,mig1d_median);

%% add transparent grey shade to the plot
for i = 0:2:n_param_update     % add shade to every other parameter
    hShade = patch([1+2*i 1+2*i, 3+2*i 3+2*i],[min(ylim) max(ylim) max(ylim) min(ylim)]...
        ,[0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha', .3);   % set FaceAlpha=.3 to make the grey shade transparent
    uistack(hShade,'bottom')    % move the shade layer to bottom 
end

%% finely tune the graph
set(gca, 'FontSize', 12)
set(gca, 'FontWeight', 'bold')
set(gca, 'XTickLabelRotation', 45)

% adjust the appearance of YLabel
hY = get(gca, 'YLabel');
hY.String = sprintf('log_1_0  parameter value');
hY.FontSize = 16;
hY.FontWeight = 'bold';

% adjust XTick and XTickLabel
nParam = length(catnames);
set(gca, 'XTick', 2*[1:nParam])
xTickLabel = catnames;

% add custom legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'Marker','d','Color',wtC);
h(2) = plot(NaN,NaN,'Marker','d','Color',gal80dC);
h(3) = plot(NaN,NaN,'Marker','d','Color',mig1dC);
legend(h, cluster_names, 'FontSize',fontsize, 'FontWeight','bold');

export_fig('violinPlotMar12th')