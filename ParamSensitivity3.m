function ParamSensitivity3(param_set, param_name, outpath)

perturb_coefficient = logspace(log10(0.01),log10(100),11);

%% load wildtype trait and a set of parameters
load('../metaData/trait_extraction/wildtype_1c.mat')
wt_trait_1c = trait;
load('../metaData/trait_extraction/wildtype_1r.mat')
wt_trait_1r = trait;

% the baseline value of the parameters
param = param_set;

% calculate the increased value and the decreased value
base_val = param.(param_name);
% vary parameter values based on the perturb coefficient
varied_value = perturb_coefficient .* base_val;

markersize = 6;

obj = struct();

%% simulate the induced level by taking the parameters into the ode
figure
set(gcf, 'position', [510 135 1122 743])

% simulate for one column fitting
sim_induce = nan(8,length(perturb_coefficient));
sim_basal = nan(8,length(perturb_coefficient));
for i = 1:length(perturb_coefficient)
    param.(param_name) = varied_value(i);
    [eval_tab, obj_value] = sensitivity_helper(param, wt_trait_1c, 'one_column');
    sim_induce(:,i) = eval_tab{:,'sim_induce'};
    sim_basal(:,i) = eval_tab{:,'sim_basal'};
    obj.col(i) = obj_value;
end
n_condition = height(eval_tab);
mid_point = ceil(length(perturb_coefficient)/2);
subplot(3,1,1)
plot(1:n_condition,sim_induce(:,1:mid_point-1), 'v', 'MarkerFaceColor', 'auto', 'markersize', markersize)
hold on
plot(1:n_condition,sim_induce(:,mid_point), 'ko', 'MarkerFaceColor', 'auto', 'markersize', markersize+3)
hold on
plot(1:n_condition,sim_induce(:,mid_point+1:end), '^', 'MarkerFaceColor', 'auto', 'markersize', markersize)

% add basal level simulation plot
hold on
plot(1:n_condition,sim_basal(:,1:mid_point-1), '<', 'markersize', markersize)
hold on
plot(1:n_condition,sim_basal(:,mid_point), 'k-', 'LineWidth', 1.5)
hold on
plot(1:n_condition,sim_basal(:,mid_point+1:end), '>', 'markersize', markersize)
% end of basal level simulation plot

hold off
set(gca, 'yscale', 'log')
grid on
xlim([0 n_condition+1])
xlabel('glucose gradient')
ylabel('sim level')
title(sprintf('%s\n', param_name))
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')

% simulate for one row fitting
sim_induce = nan(12,length(perturb_coefficient));
sim_basal = nan(12,length(perturb_coefficient));
for i = 1:length(perturb_coefficient)
    param.(param_name) = varied_value(i);
    [eval_tab, obj_value] = sensitivity_helper(param, wt_trait_1r, 'one_row');
    sim_induce(:,i) = eval_tab{:,'sim_induce'};
    sim_basal(:,i) = eval_tab{:,'sim_basal'};
    obj.row(i) = obj_value;
end
n_condition = height(eval_tab);
subplot(3,1,2)
plot(1:n_condition,sim_induce(:,1:mid_point-1), 'v', 'MarkerFaceColor', 'auto', 'markersize', markersize)
hold on
plot(1:n_condition,sim_induce(:,mid_point), 'ko', 'MarkerFaceColor', 'auto', 'markersize', markersize+3)
hold on
plot(1:n_condition,sim_induce(:,mid_point+1:end), '^', 'MarkerFaceColor', 'auto', 'markersize', markersize)

% add basal level simulation plot
hold on
plot(1:n_condition,sim_basal(:,1:mid_point-1), '<', 'markersize', markersize)
hold on
plot(1:n_condition,sim_basal(:,mid_point), 'k-', 'LineWidth', 1.5)
hold on
plot(1:n_condition,sim_basal(:,mid_point+1:end), '>', 'markersize', markersize)
% end of basal level simulation plot

hold off
set(gca, 'yscale', 'log')
grid on
xlim([0 n_condition+1])
xlabel('galactose gradient')
ylabel('sim level')
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')

%% plot obj(the sum of the deviation) versus perturb coefficient
subplot(3,1,3)
plot(1:length(perturb_coefficient), obj.col, 'LineWidth', 1.5)
hold on
plot(1:length(perturb_coefficient), obj.row, 'LineWidth', 1.5)
hold off
catnames = cell(1,length(varied_value));
for i = 1:length(varied_value)
    catnames{i} = num2str(varied_value(i), '%.4f');
end
set(gca, 'xtick', 1:length(perturb_coefficient), 'xticklabels', catnames)
grid on
legend('col fit', 'row fit', 'Location', 'best')
xlabel(sprintf('%s param values', param_name))
ylabel('sum obj')
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica')

%% export figure
export_fig(fullfile(outpath, sprintf('%s vary in logspace', param_name)));

end

function [eval_tab, obj] = sensitivity_helper( param, trait, fit_type )
output = evalGalPathway( param, trait, fit_type);

fit_type_config;

sugar_ratio = trait{index_list, 'galc'} ./ trait{index_list, 'gluc'};
eval_tab = table(output.experiment_result_linear(:,1)...
    , output.experiment_result_linear(:,2)...
    , output.simulation_result_linear(:,1)...
    , output.simulation_result_linear(:,2)...
    , trait{index_list, 'mask_basal'}...
    , trait{index_list, 'mask_induction'}...
    , sugar_ratio...
    , 'VariableNames', {'exp_basal', 'exp_induce', 'sim_basal', 'sim_induce', 'mask_basal', 'mask_induce', 'sugar_ratio'});

eval_tab = sortrows(eval_tab, 'sugar_ratio', 'ascend');
obj = output.sum_obj;

end