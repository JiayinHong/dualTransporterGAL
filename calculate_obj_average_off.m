function sum_obj = calculate_obj_average_off(trait, basal_level, induction_level)

% updated by JH on 2018.2.28, to fix the bug that when there's no OFF peak,
% sum_obj will be NaN, such as the scenario of GAL80d strain

load_global;
autofluorescence = 0;
ind_off = [trait.mask_basal] ~= 0;
ind_on = [trait.mask_induction] ~= 0;

if ind_off == 0
    obj_off = 0;
else
    % calculate expt OFF peak average and simulation OFF peak average
    expt_basal = median([trait{ind_off,'basal_level'}]);  % change mean to median to avoid extreme value
    expt_basal_linear = logyfp_to_nm(expt_basal);
    sim_basal = median([basal_level(ind_off)]);  % change mean to median to avoid extreme value
    sim_basal_linear = sim_basal + autofluorescence;
    obj_off = (log10(sim_basal_linear) - log10(expt_basal_linear))^2;
end

% calculate expt and simulation ON peak difference
expt_induce_linear = logyfp_to_nm(trait{ind_on,'ind_level'});
sim_induce_linear = induction_level(ind_on) + autofluorescence;
obj_on = (log10(sim_induce_linear) - log10(expt_induce_linear)) .^2;
obj_on = nansum(obj_on);

sum_obj = obj_off + obj_on;

end
