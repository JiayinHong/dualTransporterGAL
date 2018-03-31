function ParamSensitivity( param_set, param_name, percent_to_change )
%   this function is used to apply parameter sensitivity analysis to the
%   'param_name' from a given 'param_set', user define to what extent you
%   want to vary the parameter values, by default, it's 15%.
%   the function use the wildtype data to evaluate how the parameters fit

if nargin < 3
    fprintf('use default percent to change: 15%%\n')
    percent_to_change = 0.15;
end

% load trait
load('../metaData/trait_extraction/wildtype_1c.mat')   
wt_trait_1c = trait;
load('../metaData/trait_extraction/wildtype_1r.mat')
wt_trait_1r = trait;
load('../metaData/trait_extraction/wildtype_1r1c.mat')
wt_trait_1r1c = trait;

% the baseline value of the parameters
param = param_set;
% param_names = fieldnames(param);
figure
set(gcf, 'Position', [749 155 890 708]);
% make subplot for baseline values
subplot(3,3,1)
param_fitting_plot(param, wt_trait_1c, 'one_column');
ylabel('base value', 'FontSize', 16, 'FontWeight', 'bold');
subplot(3,3,2)
param_fitting_plot(param, wt_trait_1r, 'one_row');
subplot(3,3,3)
param_fitting_plot(param, wt_trait_1r1c, 'one_cross');

% choose the parameter you want to change, one at a time
base_val = param.(param_name);
% increase the baseline value by 'percent to change'
add_val = base_val * (1+percent_to_change);
% decrease the baseline value by 'percent to change'
min_val = base_val * (1-percent_to_change);

% make subplot for increased values
param.(param_name) = add_val;
subplot(3,3,4)
param_fitting_plot(param, wt_trait_1c, 'one_column');
ylabel(sprintf('increase by %.0f%%', 100 * percent_to_change), ...
       'FontSize', 16, 'FontWeight', 'bold');
subplot(3,3,5)
param_fitting_plot(param, wt_trait_1r, 'one_row');
subplot(3,3,6)
param_fitting_plot(param, wt_trait_1r1c, 'one_cross');

% make subplot for decreased values
param.(param_name) = min_val;
subplot(3,3,7)
param_fitting_plot(param, wt_trait_1c, 'one_column');
ylabel(sprintf('decrease by %.0f%%', 100 * percent_to_change), ...
       'FontSize', 16, 'FontWeight', 'bold');
subplot(3,3,8)
param_fitting_plot(param, wt_trait_1r, 'one_row');
subplot(3,3,9)
param_fitting_plot(param, wt_trait_1r1c, 'one_cross');

% suplabel('base line, increase, decrease', 'y');
[~,h1] = suplabel(sprintf('%s, base value = %.3e', param_name, base_val), 't');
set(h1, 'FontSize', 18)
export_fig(fullfile('../results/param_sensitivity_analysis/', sprintf('%s sensitivity analysis', param_name)));

end

