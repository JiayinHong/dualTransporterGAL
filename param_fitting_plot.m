function  [ all_conc_Glu, all_conc_Gal ] = param_fitting_plot( param, trait, fit_type )
% this version works for single gradient, one row, one column, 
% or one cross of double gradient data

output = evalGalPathway( param, trait, fit_type );
all_conc_Glu = output.all_conc_Glu;   
all_conc_Gal = output.all_conc_Gal;

% it was a bug since the sequence of output.all_conc_Glu is in the original
% one as the trait (unsorted), I fixed it by sorting all_conc_Glu in the
% sequence as sugar_ratio ascending

fit_type_config;

% switch fit_type
%     case 'one_row'
%         index_list = [4:8:92];
%     case 'one_column'
%         index_list = [65:72];
%     case 'one_cross'
%         index_list = [4:8:60,65:72,76,84,92];
%     case 'single_gradient'
%         index_list = [1:12];
% end

sugar_ratio = trait{index_list, 'galc'} ./ trait{index_list, 'gluc'};
eval_tab = table(output.experiment_result_linear(:,1)...
    , output.experiment_result_linear(:,2)...
    , output.simulation_result_linear(:,1)...
    , output.simulation_result_linear(:,2)...
    , trait{index_list, 'mask_basal'}...
    , trait{index_list, 'mask_induction'}...
    , sugar_ratio...
    , 'VariableNames', {'exp_basal', 'exp_induce', 'sim_basal', 'sim_induce', 'mask_basal', 'mask_induce', 'sugar_ratio'});

%%%%%%%% an example showing how to use sort id %%%%%%%%

% example:  [~,id] = sort(sugar_ratio);
% sorted_sugar_ratio = sugar_ratio(id);
% [~,id2] = sort(id);    % id2 is used to sort back to the original sequence
% sort_back_sugar_ratio = sorted_sugar_ratio(id2);
% 'sort_back_sugar_ratio' is exact 'sugar_ratio'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,id] = sort(eval_tab.sugar_ratio, 'ascend');
% all_conc_Glu = output.all_conc_Glu(id, :);
% all_conc_Gal = output.all_conc_Gal(id, :);

% eval_tab = sortrows(eval_tab, 'sugar_ratio', 'ascend');

n_condition = height(eval_tab);
markersize = 5;

for i_condition = 1:n_condition
    if eval_tab{i_condition, 'mask_basal'}
        plot(i_condition, eval_tab{i_condition, 'exp_basal'}, 'ok', 'markersize', markersize)
    else
        plot(i_condition, eval_tab{i_condition, 'exp_basal'}, '+k', 'markersize', markersize)
    end
    hold on
    if eval_tab{i_condition, 'mask_induce'}
        plot(i_condition, eval_tab{i_condition, 'exp_induce'}, 'or', 'markersize', markersize)
    else
        plot(i_condition, eval_tab{i_condition, 'exp_induce'}, '+r', 'markersize', markersize)
    end
end

plot(eval_tab{:,'sim_basal'}, 'k-', 'linewidth', 2)
plot(eval_tab{:,'sim_induce'}, 'r-', 'linewidth', 2)
set(gca, 'yscale', 'log')
xlim([0, n_condition+1])
title(sprintf('%s - obj : %1.4f', changeunderscore(fit_type), output.G1obj))
% title(num2str(output.sum_obj, 'obj: %1.2f'))


end

