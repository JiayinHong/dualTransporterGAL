function  [ all_conc_Glu, all_conc_Gal ] = param_fitting_plot_GAL234( param, trait, G2level, G3level, G4level, fit_type, FLAG )
% this function is modified from 'param_fitting_plot'
% this version works for plot showing how GAL3, GAL4 level fit, but not
% GAL2, since the 'eval_tab' doesn't contain G2 profiles
% the function can run normally, but not in a good structure.
% 2017.07.11 by JH

output = evalGalPathway_GAL34_changedR( param, trait, G2level, G3level, G4level, fit_type );
all_conc_Glu = output.all_conc_Glu;
all_conc_Gal = output.all_conc_Gal;

% it was a bug since the sequence of output.all_conc_Glu is in the original
% one as the trait (unsorted), I fixed it by sorting all_conc_Glu in the
% sequence as sugar_ratio ascending

load_global;
fit_type_config;
G3level = G3level(index_list,:);
G3level = logyfp_to_nm(G3level);
G4level = G4level(index_list,:);
G4level = logyfp_to_nm(G4level);

sugar_ratio = trait{index_list, 'galc'} ./ trait{index_list, 'gluc'};
eval_tab = table(G3level...
    , output.G3basal...
    , output.G3induc...
    , G4level...
    , output.G4basal...
    , output.G4induc...
    , sugar_ratio...
    , 'VariableNames', {'G3exp', 'G3sim_basal', 'G3sim_induce', 'G4exp', 'G4sim_basal', 'G4sim_induce', 'sugar_ratio'});

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
markersize = 6;

if FLAG == 'G3'
    for i_condition = 1:n_condition
        plot(i_condition, eval_tab{i_condition, 'G3exp'}, 'ok', 'markersize', markersize)
        hold on
    end
    
    plot(eval_tab{:,'G3sim_basal'}, 'k-', 'linewidth', 2)
    plot(eval_tab{:,'G3sim_induce'}, 'r-', 'linewidth', 2)
    % set(gca, 'yscale', 'log')
    xlim([0, n_condition+1])
    title(sprintf('G3 %s - obj : %1.4f', changeunderscore(fit_type), output.GAL3_obj))
end

if FLAG == 'G4'
    for i_condition = 1:n_condition
        plot(i_condition, eval_tab{i_condition, 'G4exp'}, 'ok', 'markersize', markersize)
        hold on
    end
    
    plot(eval_tab{:, 'G4sim_basal'}, 'k-', 'linewidth', 2)
    plot(eval_tab{:, 'G4sim_induce'}, 'r-', 'linewidth', 2)
    xlim([0 n_condition+1])
    title(sprintf('G4 %s - obj : %1.4f', changeunderscore(fit_type), output.GAL4_obj))
end

end

