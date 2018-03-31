function output = evalGalPathway( param, trait, fit_type, solver )
%   the third version of eval_param, aim to increase the compatibility
%   when dealing with multiple rows and columns of data
fit_type_config;
trait = trait(index_list,:);

if nargin < 4   % to be compatible with previous version
    solver = 'ode15s';  % if not specify, use matlab ode15s
end

switch solver
    case 'ode15s'
        output = getInitParamsGalPathway(param);
        y_ss_Glu = evalMultiSugarConcentrations( param, output.y0_Glu, trait.gluc, trait.galc );
        y_ss_Gal = evalMultiSugarConcentrations( param, output.y0_Gal, trait.gluc(end:-1:1), trait.galc(end:-1:1) );
        
    case 'fortran'
        output = FgetInitParamsGalPathway(param);
        y_ss_Glu = FevalMultiSugarConcentrations( param, output.y0_Glu, trait.gluc, trait.galc );
        y_ss_Gal = FevalMultiSugarConcentrations( param, output.y0_Gal, trait.gluc(end:-1:1), trait.galc(end:-1:1) );
        
end

y_ss_Gal = y_ss_Gal(end:-1:1,:);

% y_ss_Glu_init = evalMultiSugarConcentrations( param, output.y0_Glu, trait.gluc(1), trait.galc(1) );
% y_ss_Gal_init = evalMultiSugarConcentrations( param, output.y0_Gal, trait.gluc(end), trait.galc(end) );
% y_ss_Gal_init(1) = y_ss_Glu_init(1);
%
% y_ss_Glu = evalMultiSugarConcentrations( param, y_ss_Glu_init, trait.gluc(2:end), trait.galc(2:end) );
% y_ss_Glu = [y_ss_Glu_init; y_ss_Glu];
% y_ss_Gal = evalMultiSugarConcentrations( param, y_ss_Gal_init, trait.gluc(end-1:-1:1), trait.galc(end-1:-1:1) );
% y_ss_Gal = [y_ss_Gal_init; y_ss_Gal(end:-1:1,:)];

basal_level = y_ss_Glu(:,1);
induction_level = y_ss_Gal(:,1);

% [output.G1obj, output.simulation_result_linear, output.experiment_result_linear] = calculate_obj( trait, basal_level, induction_level );
[output.G1obj, ~, ~] = calculate_obj( trait, basal_level, induction_level );
% output.G1obj = calculate_obj_average_off(trait, basal_level, induction_level);

autofluorescence = get_auto_fluorescence(trait);
% output.all_conc_Glu = y_ss_Glu + autofluorescence;
% output.all_conc_Gal = y_ss_Gal + autofluorescence;
output.all_conc_Glu = y_ss_Glu; % all 12 variables concentration at steady state, initial from Glu only condition
output.all_conc_Gal = y_ss_Gal; % all 12 variables concentration at steady state, initial from Gal only condition

% for GAL1 level, add autofluorescence
output.all_conc_Glu(:,1) = output.all_conc_Glu(:,1) + autofluorescence;
output.all_conc_Gal(:,1) = output.all_conc_Gal(:,1) + autofluorescence;

end

