%% write down baseline parameter values into a file

fid = fopen(fullfile('../results/param_sensitivity_analysis/', ...
                    'param info.txt'), 'a');   
% open or create new file for writing. Append data to the end of the file.

fprintf(fid, '\n\nAppendix - baseline parameter values:\n\n');
for fdn = fieldnames(param_map)'
    fdn = fdn{1};
    fprintf(fid, '%s = %.04f\n', fdn, param.(fdn));
end

fclose(fid);

%% call function 'ParamSensitivity2' to analyze the parameters,
% vary parameter values to a larger extent
param_map = mcmc_result{1,'param_map'};
percent_to_change = 0.5;
param_obj = struct();
param_names = fieldnames(param_map);
for param_name = param_names'
    param_name = param_name{1};
    param_obj.(param_name) = ParamSensitivity2(param_map, param_name, percent_to_change);
end

%% plot the obj metric in terms of varing one parameter at a time
catnames = param_names';

n_param = length(param_names);
inc_value_cross = nan(n_param,1);   % increase param value, fit one cross data
dec_value_cross = nan(n_param,1);   % decrease param value, fit one cross data    

for i_param = 1:n_param
    inc_value_cross(i_param) = param_obj.(param_names{i_param}).cross(2);
    dec_value_cross(i_param) = param_obj.(param_names{i_param}).cross(3);
end

outpath = fullfile('../results/param_sensitivity_analysis2/', sprintf('vary_by_%.2f/', percent_to_change));
if ~exist(outpath)
    mkdir(outpath)
end

% plot obj for one cross
figure
set(gcf, 'position', [471 411 1248 459])
baseline = param_obj.(param_names{1}).cross(1);

scatter(1:n_param, inc_value_cross, 'r^', 'filled')
hold on
scatter(1:n_param, dec_value_cross, 'bv', 'filled')
hold on
plot([1 n_param], [baseline baseline], 'k-', 'LineWidth', 2)
ylim = get(gca, 'ylim');
set(gca, 'ylim', [-0.5, ylim(2)])
ylabel('obj')
set(gca, 'xtick', 1:n_param, 'xticklabels', catnames)
set(gca, 'XTickLabelRotation', 45)
grid on
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')
legend('increase by 50%','decrease by 50%')
title(sprintf('obj based on one cross fitting, percent to change = %.2f', percent_to_change))
export_fig(fullfile(outpath, sprintf('obj_one_cross, vary by %.2f', percent_to_change)));
