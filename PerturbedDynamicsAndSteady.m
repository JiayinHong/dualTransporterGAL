function PerturbedDynamicsAndSteady( origin_param, perturbed_param, trait, fit_type )
% this function is called by 'Analyze_param_sensitivity'
% the 12 species included in our model
species = {'G1', 'G2', 'G3', 'G4', 'G80', 'G3*', 'R', 'R*', 'C83', 'C84', 'glu', 'gal'};

%% plot best fitting parameter values dynamics
% solve the ode first
load_global;
opt = [];
% generate seed param for starting from GLUCOSE only state
origin_param.exglu = 0.25*perc_to_nm;
origin_param.exgal = 0;
odefunc = @(t,y)GALode(t,y,origin_param);
% when exgal==0 -> Gal3*=0, C83=0, gal=0
tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(6) = 0;     % Gal3*
tmp(9) = 0;     % C83
tmp(12) = 0;    % intracellular galactose
[t, y] = ode15s( odefunc, [0 3000], tmp, opt);

% then make the plot
figure
set(gcf, 'Position', [197 112 994 552]);
% plot dynamics starting from pure glucose
for i_species = 1:12
    subplot(3,4,i_species)
    plot(t,y(:,i_species), 'k--', 'LineWidth', 1.5)
    title(species{i_species}, 'FontSize', 15)
    hold on
end

% solve the ode initiated from galactose only state
% generate seed param for starting from GALACTOSE only state
origin_param.exglu = 0;
origin_param.exgal = 0.25*perc_to_nm;
odefunc = @(t,y)GALode(t,y,origin_param);
% when exglu==0 -> R*=0, glu=0
tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(8) = 0;     % R*
tmp(11) = 0;    % intracellular glucose
[t, y] = ode15s( odefunc, [0 3000], tmp, opt);

% plot dynamics starting from pure galactose
for i_species = 1:12
    subplot(3,4,i_species)
    plot(t,y(:,i_species), 'r--', 'LineWidth', 1.5)
    title(species{i_species}, 'FontSize', 15)
    hold on
end

legend('best fitting basal', 'best fitting induced' ...
    , 'Location', 'best')

%% plot varied parameter values dynamics

% generate seed param for starting from GLUCOSE only state
perturbed_param.exglu = 0.25*perc_to_nm;
perturbed_param.exgal = 0;
odefunc = @(t,y)GALode(t,y,perturbed_param);
% when exgal==0 -> Gal3*=0, C83=0, gal=0
tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(6) = 0;     % Gal3*
tmp(9) = 0;     % C83
tmp(12) = 0;    % intracellular galactose
[t, y] = ode15s( odefunc, [0 3000], tmp, opt);

set(gcf, 'Position', [680 234 1044 744]);
% plot dynamics starting from pure glucose
for i_species = 1:12
    subplot(3,4,i_species)
    plot(t,y(:,i_species), 'b:', 'LineWidth', 2)
    title(species{i_species}, 'FontSize', 15)
    hold on
end

% generate seed param for starting from GALACTOSE only state
perturbed_param.exglu = 0;
perturbed_param.exgal = 0.25*perc_to_nm;
odefunc = @(t,y)GALode(t,y,perturbed_param);
% when exglu==0 -> R*=0, glu=0
tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(8) = 0;     % R*
tmp(11) = 0;    % intracellular glucose
[t, y] = ode15s( odefunc, [0 3000], tmp, opt);

% plot dynamics starting from pure galactose
for i_species = 1:12
    subplot(3,4,i_species)
    plot(t,y(:,i_species), 'm:', 'LineWidth', 2)
    title(species{i_species}, 'FontSize', 15)
    hold off
end

legend('best fitting basal', 'best fitting induced' ...
    , 'perturbed basal', 'perturbed induced' ...
    , 'Location', 'best')

%% plot best fitting parameter values steady state

output = evalGalPathway( origin_param, trait, fit_type );
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

[~,id] = sort(eval_tab.sugar_ratio, 'ascend');
all_conc_Glu = output.all_conc_Glu(id, :);
all_conc_Gal = output.all_conc_Gal(id, :);
n_condition = height(eval_tab);

figure
set(gcf, 'position', [298 107 1259 822])
for i_species = 1:12
    subplot(4, 4, i_species)
    scatter(1:n_condition, all_conc_Gal(:,i_species), 'ro')
    hold on
    scatter(1:n_condition, all_conc_Glu(:,i_species), 'ko')
    hold on
    xlim([0 n_condition+1])
    set(gca, 'yscale', 'log')
    title(species{i_species}, 'FontSize', 15)
    grid on
end

% plot complex steady state
subplot(4,4,13) % G3+G3*+C83
scatter(1:n_condition, (all_conc_Gal(:,3)+all_conc_Gal(:,6)+all_conc_Gal(:,9)), 'ro')
hold on
scatter(1:n_condition, (all_conc_Glu(:,3)+all_conc_Glu(:,6)+all_conc_Glu(:,9)), 'ko')
title('G3+G3*+C83', 'FontSize', 15)
grid on
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
hold on

subplot(4,4,14) % G80+C83+C84
scatter(1:n_condition, (all_conc_Gal(:,5)+all_conc_Gal(:,9)+all_conc_Gal(:,10)), 'ro')
hold on
scatter(1:n_condition, (all_conc_Glu(:,5)+all_conc_Glu(:,9)+all_conc_Glu(:,10)), 'ko')
title('G80+C83+C84', 'FontSize', 15)
grid on
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
hold on

subplot(4,4,15) % G4+C84
scatter(1:n_condition, (all_conc_Gal(:,4)+all_conc_Gal(:,10)), 'ro')
hold on
scatter(1:n_condition, (all_conc_Glu(:,4)+all_conc_Glu(:,10)), 'ko')
title('G4+C84', 'FontSize', 15)
grid on
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
hold on

subplot(4,4,16) % R+R*
scatter(1:n_condition, (all_conc_Gal(:,7)+all_conc_Gal(:,8)), 'ro')
hold on
scatter(1:n_condition, (all_conc_Glu(:,7)+all_conc_Glu(:,8)), 'ko')
title('R+R*', 'FontSize', 15)
grid on
xlim([0 n_condition+1])
set(gca, 'yscale', 'log')
hold on

%% plot for varied parameter steady state

output = evalGalPathway( perturbed_param, trait, fit_type );
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

[~,id] = sort(eval_tab.sugar_ratio, 'ascend');
all_conc_Glu = output.all_conc_Glu(id, :);
all_conc_Gal = output.all_conc_Gal(id, :);

for i_species = 1:12
    subplot(4, 4, i_species)
    scatter(1:n_condition, all_conc_Gal(:,i_species), 'mp', 'filled')
    hold on
    scatter(1:n_condition, all_conc_Glu(:,i_species), 'bp', 'filled')
    hold on
end

% plot complex steady state
subplot(4,4,13) % G3+G3*+C83
scatter(1:n_condition, (all_conc_Gal(:,3)+all_conc_Gal(:,6)+all_conc_Gal(:,9)), 'mp')
hold on
scatter(1:n_condition, (all_conc_Glu(:,3)+all_conc_Glu(:,6)+all_conc_Glu(:,9)), 'bp')
hold off

subplot(4,4,14) % G80+C83+C84
scatter(1:n_condition, (all_conc_Gal(:,5)+all_conc_Gal(:,9)+all_conc_Gal(:,10)), 'mp')
hold on
scatter(1:n_condition, (all_conc_Glu(:,5)+all_conc_Glu(:,9)+all_conc_Glu(:,10)), 'bp')
hold off

subplot(4,4,15) % G4+C84
scatter(1:n_condition, (all_conc_Gal(:,4)+all_conc_Gal(:,10)), 'mp')
hold on
scatter(1:n_condition, (all_conc_Glu(:,4)+all_conc_Glu(:,10)), 'bp')
hold off

subplot(4,4,16) % R+R*
scatter(1:n_condition, (all_conc_Gal(:,7)+all_conc_Gal(:,8)), 'mp')
hold on
scatter(1:n_condition, (all_conc_Glu(:,7)+all_conc_Glu(:,8)), 'bp')
hold off

legend('best fitting induced', 'best fitting basal' ...
    , 'perturbed induced', 'perturbed basal' ...
    , 'Location', 'best')

end