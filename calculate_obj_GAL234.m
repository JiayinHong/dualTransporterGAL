function sum_obj = calculate_obj_GAL234( expt_level, basal_level, induction_level, basal_frac )
% this function takes in the basal level and induction level of the simulation results,
% and the experimental data (trait), then calculate the difference between them,
% return the sum of the difference as sum_obj (the objective function)

% updated by JH, 2017.07.01, based on the previous version, modified to fit
% GAL2, GAL3 & GAL4 data

load_global;    % the same formula to convert fluorescent intensity to protein level as applied to GAL1
% autofluorescence = get_auto_fluorescence(trait);
autofluorescence = 0;

% expt_basal = trait.basal_level;
% expt_induction = trait.ind_level;
experiment_result = expt_level;   % unimodal, expt data was stored in trait.basal_level

% major change as below:
simulation_result = basal_level .* basal_frac + induction_level .* (1-basal_frac);

simulation_result_linear = autofluorescence + simulation_result;

% convert the expt result from natural log to linear
experiment_result_linear = logyfp_to_nm(experiment_result);

simulation_result = log10(simulation_result_linear);
experiment_result = log10(experiment_result_linear);

sum_obj = (simulation_result - experiment_result).^2;

sum_obj = sum(sum_obj(:));

end
