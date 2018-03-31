function simulateGALPathway( trait, param_init, parameter_update, fit_type, varargin )

mcmc_version = '3.0';

% input parser
p = inputParser;

addRequired(p,'trait',@isstruct);
addRequired(p,'param_init',@isstruct);
addRequired(p,'parameter_update',@isstruct);
addRequired(p,'fit_type',@isstr);

addOptional(p,'n_propose',10000,@isnumeric);
addOptional(p,'error_tol',0.15,@isnumeric);
addOptional(p,'outfilepath','',@isstr);
addOptional(p,'jobtag','',@isstr);
addOptional(p,'arrayid','',@isstr);
addOptional(p,'allmovefrac',0,@isnumeric);

parse(p, trait, param_init, parameter_update, fit_type, varargin{:});
p = p.Results;
trait = p.trait;
param_init = p.param_init;
parameter_update = p.parameter_update;
fit_type = p.fit_type;
n_propose = p.n_propose;
error_tol = p.error_tol;
outfilepath = p.outfilepath;
jobtag = p.jobtag;
arrayid = p.arrayid;
allmovefrac = p.allmovefrac;

% setup output file and make sure path exist
if isempty(outfilepath)
    % check to see if the folder exists
%     if ~isdir('../results/mcmc_multiple_trait/')
%         mkdir('../results/mcmc_multiple_trait/')
%     end
    % Now set the file path
    task_id = str2double(arrayid);
    outfilepath = fullfile(...
        sprintf(...
        '%s-%s-%s.mat', ...
        jobtag, num2str(task_id, '%03d'), ...
        datestr(now, 'yymmdd_hhMM') ...
        ) ...
        );
end

% These next 3 lines might be deleted.
prob_data_over_parameter = struct(); % likelihood
prob_parameter_over_prior = struct(); % prior probability
param_prob = struct();              % posterior probability

for plate_names = fieldnames(trait)'
    plate_name = plate_names{1};
    % get the current trait, initial parameters, and parameter update
    trait_plate = trait.(plate_name);
    param_plate = param_init.(plate_name);
    param_update_plate = parameter_update.(plate_name);
    % accessory parameters
    n_propose_write = 5;  % the interval of iterations to output results
    
    % define function to get prior probability of a given trait and param
    get_prob_parameter_over_prior = @(param, parameter_update) get_param_prob(param, parameter_update);
    
    % init the loop
    [prob_data_over_parameter.(plate_name), simulation_result_linear] = get_prob_data_over_parameter(param_plate, trait_plate, error_tol, fit_type);
    prob_parameter_over_prior.(plate_name) = get_prob_parameter_over_prior(param_plate, param_update_plate);
    param_prob.(plate_name) = prob_data_over_parameter.(plate_name) + prob_parameter_over_prior.(plate_name);
end

sum_likelihood = prob_data_over_parameter.wt + prob_data_over_parameter.mig1d + prob_data_over_parameter.gal80d;
sum_prior_prob = prob_parameter_over_prior.wt + prob_parameter_over_prior.mig1d + prob_parameter_over_prior.gal80d;
sum_posterior_prob = param_prob.wt + param_prob.mig1d + param_prob.gal80d;

i_parameter_to_change = [];
param = param_init.wt;

% set up variables to record results
accept_list = nan(n_propose, 1);
sum_likelihood_list = nan(n_propose, 1);
sum_prior_prob_list = nan(n_propose, 1);
sum_posterior_prob_list = nan(n_propose, 1);
running_t_list = nan(n_propose, 1);
param_list = struct();
for fdn = fieldnames(param)'
    param_list.(fdn{1}) = nan(n_propose, 1);
end
i_param_list = 1;
param_map = param;  % map estimate
param_prob_map = sum_posterior_prob; % map posterior probability
parameter_updated_list = cell(n_propose, 1);
simulation_result_linear_list = cell(n_propose, 1);

% prep the output
save(outfilepath ...
    , 'trait' ...
    , 'param_init' ...
    , 'parameter_update' ...
    , 'error_tol' ...
    , 'outfilepath' ...
    , 'jobtag' ...
    , 'n_propose' ...
    , 'mcmc_version' ...
    );

% MAIN LOOP
tic
n_parameter_update = height(parameter_update.wt);
for i_propose = 1:n_propose
    if rand < allmovefrac
        % block move, special case, all paramters moves
        flag_block_move = 1;
        i_parameter_update_list = 1:n_parameter_update;
        parameter_updated_label = 'all';
    else
        % single parameter move
        flag_block_move = 0;
        if isempty(i_parameter_to_change)
            % the beginning of a new update cycle
            
            i_parameter_to_change = randperm(n_parameter_update);
            % generate a random permutation of the parameters to be updated
        end
        i_parameter_update_list = i_parameter_to_change(end);
        % each time pop out the last parameter of the permutation and
        % change it
        i_parameter_to_change(end) = []; % then remove the last param from the permutation
        parameter_updated_label = parameter_update.wt{i_parameter_update_list(1), 'parameter_name'}{1};
    end
    
    % create new parameter
    param_new = struct;
    param_new.wt = param;
    for i_parameter_update = i_parameter_update_list
        parameter_update_name = parameter_update.wt{i_parameter_update, 'parameter_name'}{1};
        if rand < parameter_update.wt{i_parameter_update, 'gaussian_mix_fraction'}
            proposal_scale = parameter_update.wt{i_parameter_update, 'proposal_sigma_1'};
        else
            proposal_scale = parameter_update.wt{i_parameter_update, 'proposal_sigma_2'};
        end
        param_new.wt.(parameter_update_name) = random(...
            'lognormal' ...
            , log(param_map.(parameter_update_name))...
            , proposal_scale ...
            );
    end
    
    param_new.mig1d = param_new.wt;
    param_new.mig1d.aR = 0;
    
    param_new.gal80d = param_new.wt;
    param_new.gal80d.a80 = 0;
    param_new.gal80d.ag80 = 0;
    
    % remember value of moved paramter
    if flag_block_move
        old_val = nan; % block move
        new_val = nan;
    else
        old_val = param.(parameter_update_name);
        new_val = param_new.wt.(parameter_update_name);
    end
    
    % evaluate the prob
    for plate_names = fieldnames(trait)'
        plate_name = plate_names{1};
        [prob_data_over_parameter_new.(plate_name), simulation_result_linear_new] = get_prob_data_over_parameter(param_new.(plate_name), trait.(plate_name), error_tol, fit_type);
        prob_parameter_over_prior_new.(plate_name) = get_prob_parameter_over_prior(param_new.(plate_name), parameter_update.(plate_name));
        param_prob_new.(plate_name) = prob_data_over_parameter_new.(plate_name) + prob_parameter_over_prior_new.(plate_name);
    end
    
    sum_likelihood_new = prob_data_over_parameter_new.wt + prob_data_over_parameter_new.mig1d + prob_data_over_parameter_new.gal80d;
    sum_prior_prob_new = prob_parameter_over_prior_new.wt + prob_parameter_over_prior_new.mig1d + prob_parameter_over_prior_new.gal80d;
    sum_posterior_prob_new = param_prob_new.wt + param_prob_new.mig1d + param_prob_new.gal80d;
    
    % evaluate jumping prob
    relative_transition_prob = 0;  % symmetric move in log scale, hence 0
    
    % evaluate whether accept or not
    accept_thre = sum_posterior_prob_new - sum_posterior_prob + relative_transition_prob;
    accept = random('uniform', 0, 1, 1) < exp(accept_thre);
    if accept
        param = param_new.wt;
        sum_likelihood = sum_likelihood_new;
        sum_prior_prob = sum_prior_prob_new;
        sum_posterior_prob = sum_posterior_prob_new;
        simulation_result_linear = simulation_result_linear_new;
        
        if sum_posterior_prob > param_prob_map
            param_prob_map = sum_posterior_prob;
            param_map = param;
        end
    end
    
    accept_list(i_propose) = accept;
    sum_likelihood_list(i_propose) = sum_likelihood;
    sum_prior_prob_list(i_propose) = sum_prior_prob;
    sum_posterior_prob_list(i_propose) = sum_posterior_prob;
    running_t_list(i_propose) = toc;
    parameter_updated_list{i_propose} = parameter_updated_label;
    simulation_result_linear_list{i_propose} = simulation_result_linear;
    
    % screen output
    fprintf('%d\t%d\t%2.2f\t%.2f\t%.2f\t%.2f\t%s\t%1.2e\t%1.2e\n', ...
        i_propose, ...
        accept, ...
        -sum_likelihood * error_tol^2, ...
        sum_likelihood, ...
        sum_prior_prob, ...
        sum_posterior_prob, ...
        parameter_updated_label, ...
        old_val, new_val);
    toc
    
    for fdn = fieldnames(param_list)'
        param_list.(fdn{1})(i_param_list) = param.(fdn{1});
    end
    i_param_list = i_param_list+1;
    
    if (i_propose < 100 && mod(i_propose, n_propose_write) == 0) ...
            || (i_propose >= 100 && mod(i_propose, 100*n_propose_write) == 0) ...
            || i_propose == n_propose
        
        save(outfilepath ...
            , 'param_list' ...
            , 'sum_likelihood_list' ...
            , 'sum_prior_prob_list' ...
            , 'sum_posterior_prob_list' ...
            , 'accept_list' ...
            , 'param_prob_map' ...
            , 'param_map' ...
            , 'i_propose' ...
            , 'running_t_list' ...
            , 'simulation_result_linear_list' ...
            , 'parameter_updated_list' ...
            , '-append'...
            );
    end
    
end

end


function res = get_param_prob(param, parameter_update)
res = 0;
for i = 1:height(parameter_update)
    res = res + ...
        + log(...
        pdf(...
        'lognormal', ...
        param.(parameter_update{i, 'parameter_name'}{1}), ...
        log(parameter_update{i, 'prior_mean'}), ...
        parameter_update{i, 'prior_sigma'} ...
        ) ...
        );
end
end

function [prob, simulation_result_linear] = get_prob_data_over_parameter(param, trait, error_tol, fit_type)

output = evalGalPathway(param, trait, fit_type);
prob = - output.sum_obj / error_tol^2;
simulation_result_linear = output.simulation_result_linear;

end