if 0
    %% heatmap, cluster & dendrogram for single mutants data
    
    % load obj results from each file, then store them all in 'obj_all'
    clear
    tmp = dir('../paramScanResults/');
    tmp = tmp([tmp.isdir]~=1);  % filter out the directory
    obj_all = struct;
    for i=1:length(tmp)
        if ~strcmp(tmp(i).name, '.DS_Store')
            load(fullfile('../paramScanResults/',tmp(i).name))
            fdn = fieldnames(obj);
            fdn = fdn{1};
            obj_all.(fdn) = obj.(fdn);
        end
    end
    
    % store all the strain names into a cell array
    clear strNames
    nStr = numel(fieldnames(obj_all));
    strNames = cell(nStr,1);
    
    for i=1:length(allTraits)
        tmp = allTraits(i).name;
        strNames{i} = tmp(1:end-4);   % to remove '.mat'
    end
    
    % organize the data
    nPerturb = size(obj_all.str01,1);
    mydata = nan(nStr,nPerturb);
    for iStr = 1:nStr
        fdn = num2str(sprintf('str%02d',iStr));
        for iPerturb = 1:nPerturb
            mydata(iStr,iPerturb) = nanmin(obj_all.(fdn)(iPerturb,:));
        end
    end
    
    % make a list of sorted strain names
%     load('../wildtype_96well-025-170925_16:18.mat')
%     param_names = parameter_update.parameter_name;
%     rowLabels = param_names;
%     
%     for i=1:length(param_names)
%         if regexp(param_names{i},'n*')
%             tmp1 = sprintf([param_names{i}, '=4']);
%             tmp2 = sprintf([param_names{i}, '=1']);
%         else
%             tmp1 = sprintf([param_names{i}, '*10']);
%             tmp2 = sprintf([param_names{i}, '*0.1']);
%         end
%         sortedNames{2*i-1} = tmp1;
%         sortedNames{2*i} = tmp2;
%     end
%     
%     [~, ind] = ismember(sortedNames, strNames);
%     mySortedData = mydata(ind,:);
%     colLabels = sortedNames;
    load('best_param.mat','param')
    param_names = fieldnames(param);
    rowLabels = param_names;
    colLabels = strNames;
    
    % make heatmap
    Cmap = flipud(parula);
    
    figure
    set(gcf, 'position', [245 101 1628 756])
    h=heatmap(rowLabels, colLabels, mydata, 'Colormap', Cmap ...
        , 'ColorScaling', 'log', 'CellLabelFormat', '%.0f' ...
        , 'ColorLimits', [log(7), log(100)], 'FontSize', 12);
    
    % heatmap with bi-directional cluster and dendrogram
    % CG = clustergram(mySortedData, 'RowLabels',colLabels, 'ColumnLabels',rowLabels ...
    %     , 'Standardize','none', 'Cluster','all', 'Colormap',Cmap ...
    %     , 'LogTrans', 1);
    % cgAxes = plot(CG);
    % set(gcf, 'position', [360 50 569 655])
    % % 'LogTrans' transform the data from natural scale to log2 scale
    % set(cgAxes, 'Clim', [log2(.1), log2(100)], 'FontSize', 12)
    
    %% plot obj variation by single parameter scan
    
    param_name = 'aR';
    idx = strfind(param_names,param_name);
    idx = find(~cellfun(@isempty,idx));
    obj = obj_all.str02(idx,:);

    % the baseline value of the parameters
    load('best_param.mat','param')
    base_val = param.(param_name);  
    
    if regexp(param_name, 'n*')     % hill coefficients
        varied_value = 1:13;
    else
    perturb_coefficient = logspace(log10(0.001), log10(1000), 49);
    % vary parameter values based on the perturb coefficient
    varied_value = perturb_coefficient .* base_val;
    end
    
    nValue = 49;
    obj_thresh = 10;
    
    % plot obj(the sum of the deviation) versus perturb coefficient
    figure
    set(gcf, 'position', [2049 418 1065 411])
    
    plot(1:nValue, obj, 'LineWidth', 1.5, 'Marker', 'o')
    
    catnames = cell(1,nValue);
    for i = 1:2:nValue
        catnames{i} = num2str(varied_value(i), '%.4f');
%         catnames{i} = num2str(varied_value(i), '%.3e');
    end
    set(gca, 'xtick', 1:nValue, 'xticklabels', catnames)
    set(gca, 'XTickLabelRotation', 45)
    grid on
    
    hline = refline([0 obj_thresh]);    % add a reference line
    hline.LineWidth = 1.5;
    hline.Color = 'm';
    hline.LineStyle = '-.';
    
    xlabel(sprintf('%s param values', param_name))
    ylabel('sum obj - 96 well')
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')
    
    %% plot a dendrogram tree to view the dissimilarity between the perturbations
    
    % mySortedData is Strains * Perturbations, since we'd like to know which
    % clusters of perturbations have similar effect to the system, so we need
    % to first transpose the matrix, and then calculate the pairwise distance
    % of each two perturbations using default metric - Euclidean distance
    distList = pdist(mySortedData');
    
    % the 1st and 2nd column of the linkaged matrix contains the cluster indices
    % linked in pairs to form a binary tree, while the 3rd column contains the
    % linkage distances between the two clusters merged in the same row
    tree = linkage(distList);
    
    % dendrogram is a way to visualize the linkage tree, while the x tick
    % labels are the original indices of the objects, the height of the
    % U-shaped lines indicates the distance between the objects
    figure
    set(gcf,'position',[680 511 950 467])
    H = dendrogram(tree);
    for i=1:numel(H)
        H(i).LineWidth = 1.5;
    end
    
    perturbID = str2num(get(gca,'XTickLabels'));
    set(gca, 'XTickLabels', param_names(perturbID))
    set(gca, 'XTickLabelRotation', 45, 'FontSize', 13)
    
    
    %% heatmap, cluster & dendrogram for multiple strains data
    
    % load obj results from each file, then store them all in 'obj_all'
    clear
    obj_all = struct;
    tmp = dir('../multiStrainFromKayla/');
    for i=3:length(tmp)
        if ~ strcmp(tmp(i).name, '.DS_Store')
            load(fullfile('../multiStrainFromKayla/', tmp(i).name))
            str = fieldnames(obj);
            str = str{1};
            no = str(4:end);
            no = str2num(no)-2;
            fdn = sprintf('str%02d',no);
            obj_all.(fdn) = obj.(str);
        end
    end
    
    % store all the strain names into a cell array
    clear strNames
    nStr = numel(fieldnames(obj_all));
    strNames = cell(nStr,1);
    
    for i=3:length(allTraits)
        if ~ strcmp(allTraits(i).name, '.DS_Store')
            tmp = allTraits(i).name;
            strNames{i-2} = tmp(1:end-4);
        end
    end
    
    % organize the data
    nPerturb = size(obj_all.str01,1);
    mydata = nan(nStr,nPerturb);
    for iStr = 1:nStr
        fdn = num2str(sprintf('str%02d',iStr));
        for iPerturb = 1:nPerturb
            mydata(iStr,iPerturb) = min(obj_all.(fdn)(iPerturb,:));
        end
    end
    
    % hard-coded the list of picked strains
    load('../wildtype_96well-025-170925_16:18.mat')
    param_names = parameter_update.parameter_name;
    rowLabels = param_names;
    
    pickedStrains = strNames([1,2,4:6,8,10:13,16,19,20,22,23,26,27,29,31:34],:);
    
    [~,ind] = ismember(pickedStrains, strNames);
    myPickedData = mydata(ind,:);
    colLabels = pickedStrains;
    
    % make heatmap
    Cmap = flipud(parula);
    Climits = [min(myPickedData(:)), max(myPickedData(:))];
    % Z_Climits = [min(zscore(myPickedData(:))), max(zscore(myPickedData(:)))];
    
    % heatmap without cluster or dendrogram tree
    figure
    set(gcf, 'position', [218 129 780 496])
    heatmap(rowLabels, colLabels, myPickedData, 'Colormap', Cmap ...
        , 'ColorLimits', Climits, 'FontSize', 12);
    
    % heatmap with bi-directional cluster and dendrogram
    CG = clustergram(myPickedData, 'RowLabels',colLabels, 'ColumnLabels',rowLabels ...
        , 'Standardize','none', 'Cluster','all', 'Colormap',Cmap);
    cgAxes = plot(CG);
    set(gcf, 'position', [397 99 810 552])
    set(cgAxes, 'Clim', Climits, 'FontSize', 12)
end

%% visualization of single parameter scan results
% batch simulated on Nov.28
clear

tmp = dir('../paramScanResults/');
tmp = tmp([tmp.isdir]==1);  % to get rid of '.DS_Store'
for i=1:length(tmp)
    if ~ismember(tmp(i).name,{'.','..'})
    storeDir{i} = ['../paramScanResults/',tmp(i).name,'/'];
    end
end
storeDir = storeDir(~cellfun('isempty',storeDir));

for input_path = storeDir
    input_path = input_path{1};
    clusterDendrogramHelper(input_path)
end


%% below is in-script function
function clusterDendrogramHelper(input_path)
obj_all = struct;
tmp = dir(input_path);
tmp = tmp([tmp.isdir]==0);
for i=1:length(tmp)
    if ~ strcmp(tmp(i).name, '.DS_Store')
        load(fullfile(input_path, tmp(i).name))
        str = fieldnames(obj);
        str = str{1};
        obj_all.(str) = obj.(str);
    end
end

% store all the strain names into a cell array
clear strNames
nStr = numel(fieldnames(obj_all));
strNames = cell(nStr,1);

for i=1:length(allTraits)
    if ~ strcmp(allTraits(i).name, '.DS_Store')
        tmp = allTraits(i).name;
        strNames{i} = tmp(1:end-4);
    end
end

% organize the data
nPerturb = size(obj_all.str01,1);
mydata = nan(nStr,nPerturb);
fdnAll = fieldnames(obj_all);
for iStr = 1:nStr
    fdn = fdnAll{iStr};
    for iPerturb = 1:nPerturb
        mydata(iStr,iPerturb) = min(obj_all.(fdn)(iPerturb,:));
    end
end

% setup row labels and col labels
parameter_update = readtable('Nov10th_param_config_set10.csv');
param_names = parameter_update.parameter_name;
rowLabels = param_names;
colLabels = strNames;

% make heatmap
Cmap = flipud(parula);
Climits = [min(mydata(:)), max(mydata(:))];
% Z_Climits = [min(zscore(myPickedData(:))), max(zscore(myPickedData(:)))];

tok = regexp(input_path,'.*/(\w+-\d{3})/','tokens');
tok = tok{1};

% heatmap without cluster or dendrogram tree
figure
set(gcf, 'position', [218 129 780 496])
h = heatmap(rowLabels, colLabels, mydata, 'Colormap', Cmap ...
    , 'ColorLimits', Climits, 'FontSize', 12);
h.Title = ['Starting from ''',tok{1},''''];
export_fig(fullfile(input_path,h.Title))

% heatmap with bi-directional cluster and dendrogram
CG = clustergram(mydata, 'RowLabels',changeunderscore(colLabels), 'ColumnLabels',rowLabels ...
    , 'Standardize','none', 'Cluster','all', 'Colormap',Cmap);
cgAxes = plot(CG);
set(gcf, 'position', [397 99 810 552])
set(cgAxes, 'Clim', Climits, 'FontSize', 12)
cgAxes.Title.String = ['Starting from ''',changeunderscore(tok{1}),''''];
cgAxes.Title.Position = cgAxes.Title.Position + [0 8 0];
export_fig(fullfile(input_path, [h.Title,' with dendrogram']))

end