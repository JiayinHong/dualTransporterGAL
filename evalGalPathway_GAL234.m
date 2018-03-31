function output = evalGalPathway_GAL234( param, trait, G2level, G3level, G4level, fit_type )
% this function is modified from 'evalGalPathway'
% it's called by 'mcmc_for_GAL234' and 'param_fitting_plot_GAL234'
% the function calculates the objective function value of GAL1, GAL3, GAL4
% separately, then sum up them together.
% to calculate GAL1 obj, it calls 'calculate_obj', while for the others, it
% calls 'calculate_obj_GAL234'. 
% 2017.07.11 by JH

output = getInitParamsGalPathway(param);

fit_type_config;
trait = trait(index_list,:);

y_ss_Glu = evalMultiSugarConcentrations( param, output.y0_Glu, trait.gluc, trait.galc );
y_ss_Gal = evalMultiSugarConcentrations( param, output.y0_Gal, trait.gluc(end:-1:1), trait.galc(end:-1:1) );
y_ss_Gal = y_ss_Gal(end:-1:1,:);


% first, calculate obj for GAL1 fitting
G1basal = y_ss_Glu(:,1);
G1induc = y_ss_Gal(:,1);

[GAL1_obj, ~, ~] = calculate_obj( trait, G1basal, G1induc );
basal_frac = trait.basal_frac;

% then, calculate obj for GAL2 fitting
if G2level == 0     % do not fit G2
    GAL2_obj = 0;   % July 1st tried to fit all the GAL1,2,3 data, which turned out to be very bad results,
                    % thus try if only fit GAL1 and GAL3, is it possible to get a good fitting
else
    G2level = G2level(index_list,:);
    G2basal = y_ss_Glu(:,2);
    G2induc = y_ss_Gal(:,2);
    GAL2_obj = calculate_obj_GAL234( G2level, G2basal, G2induc, basal_frac );
end

% then, calculate obj for GAL3 fitting
if G3level == 0     % do not fit G3
    output.G3basal = 0;
    output.G3induc = 0;
    GAL3_obj = 0;
else
    G3level = G3level(index_list,:);
    
    % G3basal = y_ss_Glu(:,3);
    % G3induc = y_ss_Gal(:,3);
    G3basal = y_ss_Glu(:,3) + y_ss_Glu(:,6) + y_ss_Glu(:,9);    % G3, G3*, C83
    G3induc = y_ss_Gal(:,3) + y_ss_Gal(:,6) + y_ss_Gal(:,9);
    % the experimental data is GAL3 promoter driven fluorescent proteins, which
    % should be the total amount of every species that contains a GAl3
    % component, no matter free or combined
    output.G3basal = G3basal;
    output.G3induc = G3induc;
    GAL3_obj = calculate_obj_GAL234( G3level, G3basal, G3induc, basal_frac );   
end

% then, calculate obj for GAL4 fitting
if G4level == 0     % do not fit G4
    output.G4basal = 0;
    output.G4induc = 0;
    GAL4_obj = 0;
else 
    G4level = G4level(index_list,:);
    G4basal = y_ss_Glu(:,4) + y_ss_Glu(:,10);       % G4, C84
    G4induc = y_ss_Gal(:,4) + y_ss_Gal(:,10);
    output.G4basal = G4basal;
    output.G4induc = G4induc;
    GAL4_obj = calculate_obj_GAL234( G4level, G4basal, G4induc, basal_frac );
end

% finally, sum up all three obj
output.GAL3_obj = GAL3_obj;
output.GAL4_obj = GAL4_obj;
output.sum_obj = GAL1_obj+GAL2_obj+GAL3_obj+GAL4_obj;

output.all_conc_Glu = y_ss_Glu; % all 12 variables concentration at steady state, initial from Glu only condition
output.all_conc_Gal = y_ss_Gal; % all 12 variables concentration at steady state, initial from Gal only condition

end

