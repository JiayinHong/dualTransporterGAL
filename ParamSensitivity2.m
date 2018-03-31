function obj = ParamSensitivity2( param_set, param_name, percent_to_change )
%   This function aims to compare three metrics between baseline parameter
%   vales and increased or decreased parameter values, the induced level,
%   the objective function, and the decision threshold. For the induced
%   level, each subplot show the scenario for one parameter variation, use
%   marker circle to represent the level of baseline value, and
%   upward-pointing triangle to represent the level of increased value, and
%   downward-pointing triangle to represent the level of decrease value
%   the output 'obj' has three fields: col, row, and cross. For each field,
%   the first one is the baseline value, the second one is the increased
%   value, while the third one is the decreased value.

% 2018.01.08 updated by JH
% update the script to analyze parameter sensitivity in regards of fitting
% BC187 strain, one-cross

%% load BC187 trait and a set of parameters
% load('../metaData/trait_extraction/BC187_Kayla_Nov29.mat')
% bc_trait = trait;
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat')

% the baseline value of the parameters
param = param_set;

% calculate the increased value and the decreased value
base_val = param.(param_name);
% increase the baseline value by 'percent to change'
add_val = base_val * (1+percent_to_change);
% decrease the baseline value by 'percent to change'
min_val = base_val * (1-percent_to_change);

obj = struct();

%% simulate the induced level by taking the parameters into the ode

% simulate for one cross fitting
param.(param_name) = base_val;
[base_obj] = sensitivity_helper(param, trait, '96well');
obj.cross(1) = base_obj;

param.(param_name) = add_val;
[inc_obj] = sensitivity_helper(param, trait, '96well');
obj.cross(2) = inc_obj;

param.(param_name) = min_val;
[dec_obj] = sensitivity_helper(param, trait, '96well');
obj.cross(3) = dec_obj;

end


function [obj] = sensitivity_helper( param, trait, fit_type )
output = evalGalPathway( param, trait, fit_type);

fit_type_config;

sugar_ratio = trait{index_list, 'galc'} ./ trait{index_list, 'gluc'};
% eval_tab = table(output.experiment_result_linear(:,1)...
%     , output.experiment_result_linear(:,2)...
%     , output.simulation_result_linear(:,1)...
%     , output.simulation_result_linear(:,2)...
%     , trait{index_list, 'mask_basal'}...
%     , trait{index_list, 'mask_induction'}...
%     , sugar_ratio...
%     , 'VariableNames', {'exp_basal', 'exp_induce', 'sim_basal', 'sim_induce', 'mask_basal', 'mask_induce', 'sugar_ratio'});
% 
% eval_tab = sortrows(eval_tab, 'sugar_ratio', 'ascend');
obj = output.G1obj;

end