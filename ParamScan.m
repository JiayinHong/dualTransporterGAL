function ParamScan(param_set, param_name, aim_trait, outpath)
% caution: input 'param_name' should be a cell array

% this function is used to scan single parameter values and joint parameter
% values in mig1d & gal80d strains, starting from parameter values of
% "best-fit" of wildtype, the goal is to see by solely tuning one or two
% parameters, is it even possible to fit mutants from a good fit of
% wildtype

% this function first calculate the difference in each of 96 wells, then
% directly make the plot using the values. For batch analysis, it's better
% to store the calculation results somewhere, see 'ParamScan2'

param = param_set;
perturb_coefficient = logspace(log10(0.01),log10(100),9);
nValue = length(perturb_coefficient);

%% scan joint parameter values
if numel(param_name)==2
    baseParam1 = param.(param_name{1});
    baseParam2 = param.(param_name{2});
    varyParam1 = perturb_coefficient .* baseParam1;
    varyParam2 = perturb_coefficient .* baseParam2;
    rowLabels = num2str(varyParam2(:), '%.3f');
    colLabels = num2str(varyParam1(:), '%.3f');
    
    tic
    for i=1:nValue
        param.(param_name{1}) = varyParam1(i);
        
        for j=1:nValue
            param.(param_name{2}) = varyParam2(j);
            
            output = evalGalPathway(param, aim_trait, '96well');
            obj(i,j) = output.G1obj;
        end
    end
    toc
    
    % R2016a version
    figure
    set(gcf, 'position', [689 136 1036 811])
    heatmap(obj, rowLabels, colLabels, '%.2f', 'Colorbar', true ...
        , 'Colormap', flipud(hot) ...       % to flip the colormap so that lower numerical values are brighter
        , 'MinColorValue', 1, 'MaxColorValue', 300 ...  % to normalize among different heatmaps
        , 'ShowAllTicks', true, 'TextColor', 'r', 'FontSize', 14  ...
        , 'GridLines', ':', 'ColorLevels', 128, 'TickFontSize', 15);
    h=title(sprintf('%s & %s parameter value scan', param_name{1}, param_name{2}), 'FontSize', 15);
    xlabel(param_name{2})
    ylabel(param_name{1})
    export_fig(fullfile(outpath, sprintf('%s & %s parameter value scan', param_name{1}, param_name{2})));
    
    % R2017a version
%     figure
%     set(gcf, 'position', [689 136 1036 811])
%     heatmap(rowLabels, colLabels, obj, 'CellLabelFormat', '%.2f', 'FontSize', 12);
%     title(sprintf('%s & %s parameter value scan', param_name{1}, param_name{2}));
%     export_fig(fullfile(outpath, get(gca, 'Title')))
    
    %% scan single parameter value
else
    param_name = param_name{1};
    perturb_coefficient = logspace(log10(0.01),log10(100),11);
%     perturb_coefficient = logspace(log10(0.001), log10(1000), 13);
    nValue = length(perturb_coefficient);
    obj_thresh = 10;
    
    % the baseline value of the parameters
    base_val = param.(param_name);   
    % vary parameter values based on the perturb coefficient
    varied_value = perturb_coefficient .* base_val;
    
    if regexp(param_name, 'n*')     % hill coefficients
        perturb_coefficient = linspace(1,4,7);
        nValue = length(perturb_coefficient);
        varied_value = perturb_coefficient;
    end    
    
    obj = [];
    
    for i=1:nValue
        param.(param_name) = varied_value(i);
        output = evalGalPathway(param, aim_trait, '96well');
        obj(i) = output.G1obj;
    end
    
    if min(obj)<obj_thresh
        fprintf('\nSolely tuning %s is possible to fit mutant data\n', param_name)
    end
    
    % plot obj(the sum of the deviation) versus perturb coefficient
    figure
    set(gcf, 'position', [2049 418 1065 411])
    
    plot(1:nValue, obj, 'LineWidth', 1.5, 'Marker', 'o')
    
    catnames = cell(1,nValue);
    for i = 1:nValue
        catnames{i} = num2str(varied_value(i), '%.4f');
%         catnames{i} = num2str(varied_value(i), '%.3e');
    end
    set(gca, 'xtick', 1:nValue, 'xticklabels', catnames)
%     set(gca, 'XTickLabelRotation', 45)
    grid on
    
    hline = refline([0 obj_thresh]);    % add a reference line
    hline.LineWidth = 1.5;
    hline.Color = 'm';
    hline.LineStyle = '-.';
    
    xlabel(sprintf('%s param values', param_name))
    ylabel('sum obj - 96 well')
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')
    
    % export figure
    export_fig(fullfile(outpath, sprintf('%s vary in logspace', param_name)));
    
end
end
