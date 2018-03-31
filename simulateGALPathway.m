function simulateGALPathway(trait, param_init, parameter_update, fit_type, varargin)

mcmc_version = '3.0';

% input parser
p = inputParser;

addRequired(p, 'trait', @isstruct);
addRequired(p, 'param_init', @isstruct);
addRequired(p, 'parameter_update', @istable);
addRequired(p, 'fit_type', @isstr);

addOptional(p, 'n_propose', 50000, @isnumeric);
addOptional(p, 'error_tol', 0.15, @isnumeric);
addOptional(p, 'outfilepath', '', @isstr);
addOptional(p, 'jobtag', '', @isstr);
addOptional(p, 'arrayid', '', @isstr);
addOptional(p, 'allmovefrac', 0.2, @isnumeric);

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

% setup output filename
if isempty(outfilepath)
	outfilepath = fullfile(...
		sprintf(...
		'%s-%03d-%s.mat', ...
		jobtag, str2double(arrayid), datestr(now, 'yymmdd_hhMM')...
		)...
		);
end

% accessory parameters
n_propose_write = 5;	% record results to output file every 5 iterations

param_wt = param_init;
param_gal80d = param_init;
param_gal80d.a80 = 0;
param_gal80d.ag80 = 0;
param_mig1d = param_init;
param_mig1d.aR = 0;

trait_wt = trait.wildtype;
trait_gal80d = trait.gal80d;
trait_mig1d = trait.mig1d;

% define function to get prior probability of a given parameter value
get_prob_parameter_over_prior = @(param) get_param_prob(param, parameter_update);

% init the loop
wt_prob_data_over_parameter = get_prob_data_over_parameter(param_wt, trait_wt, error_tol, fit_type);
gal80d_prob_data_over_parameter = get_prob_data_over_parameter(param_gal80d, trait_gal80d, error_tol, fit_type);
mig1d_prob_data_over_parameter = get_prob_data_over_parameter(param_mig1d, trait_mig1d, error_tol, fit_type);
sum_likelihood = wt_prob_data_over_parameter + gal80d_prob_data_over_parameter + mig1d_prob_data_over_parameter;

prior_prob = get_prob_parameter_over_prior(param_init);

posterior_prob = sum_likelihood + prior_prob;

param = param_init;
param_map = param_init;
param_prob_map = posterior_prob;

% setup variables to record results
accept_list = nan(n_propose, 1);
wt_likelihood_list = nan(n_propose, 1);
gal80d_likelihood_list = nan(n_propose, 1);
mig1d_likelihood_list = nan(n_propose, 1);
sum_likelihood_list = nan(n_propose, 1);
prior_prob_list = nan(n_propose, 1);
posterior_prob_list = nan(n_propose, 1);
running_t_list = nan(n_propose, 1);
parameter_updated_list = cell(n_propose, 1);
param_list = struct();
for fdn = fieldnames(param_init)'
	param_list.(fdn{1}) = nan(n_propose, 1);
end

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
n_parameter_update = height(parameter_update);
i_parameter_to_change = [];
i_param_list = 1;

for i_propose = 1:n_propose
	if rand < allmovefrac	% block move, special case, all parameters vary
		flag_block_move = 1;
		i_parameter_update_list = 1:n_parameter_update;
		parameter_updated_label = 'all';
	else
		flag_block_move = 0;	% single parameter vary
		if isempty(i_parameter_to_change)	% the beginning of a new update cycle
			i_parameter_to_change = randperm(n_parameter_update);	% generate a random permutation of the parameters to vary
		end
		i_parameter_update_list = i_parameter_to_change(end);
		% each time pop out the last parameter of the permutation and change it
		i_parameter_to_change(end) = [];	% then remove the last param from the permutation
		parameter_updated_label = parameter_update{i_parameter_update_list, 'parameter_name'}{1};
	end
	
	% create new parameter
	param_new = param;
	for i_parameter_update = i_parameter_update_list
		parameter_update_name = parameter_update{i_parameter_update, 'parameter_name'}{1};
		if rand < parameter_update{i_parameter_update, 'gaussian_mix_fraction'}
			proposal_scale = parameter_update{i_parameter_update, 'proposal_sigma_1'};
		else
			proposal_scale = parameter_update{i_parameter_update, 'proposal_sigma_2'};
		end
		param_new.(parameter_update_name) = random(...
			'lognormal', log(param_map.(parameter_update_name)), proposal_scale);
	end
	
	param_new_gal80d = param_new;
	param_new_gal80d.a80 = 0;
	param_new_gal80d.ag80 = 0;
	
	param_new_mig1d = param_new;
	param_new_mig1d.aR = 0;
	
	% record value of moved parameter
	if flag_block_move
		old_val = nan;
		new_val = nan;	% no screen display if block move
	else
		old_val = param.(parameter_update_name);
		new_val = param_new.(parameter_update_name);
	end
	
	% evaluate the prob
	wt_new_prob_data_over_parameter = get_prob_data_over_parameter(param_new, trait_wt, error_tol, fit_type);
	gal80d_new_prob_data_over_parameter = get_prob_data_over_parameter(param_new_gal80d, trait_gal80d, error_tol, fit_type);
	mig1d_new_prob_data_over_parameter = get_prob_data_over_parameter(param_new_mig1d, trait_mig1d, error_tol, fit_type);
	sum_likelihood_new = wt_new_prob_data_over_parameter + gal80d_new_prob_data_over_parameter + mig1d_new_prob_data_over_parameter;

	prior_prob_new = get_prob_parameter_over_prior(param_new);

	posterior_prob_new = sum_likelihood_new + prior_prob_new;
	
	% evaluate whether accept or not
	accept_thre = posterior_prob_new - posterior_prob;
	accept = random('uniform', 0, 1, 1) < exp(accept_thre);
	if accept
		param = param_new;
		sum_likelihood = sum_likelihood_new;
		prior_prob = prior_prob_new;
		posterior_prob = posterior_prob_new;
		
		if posterior_prob > param_prob_map
			param_prob_map = posterior_prob;
			param_map = param;
		end
	end
	
	accept_list(i_propose) = accept;
	wt_likelihood_list(i_propose) = wt_new_prob_data_over_parameter;
	gal80d_likelihood_list(i_propose) = gal80d_new_prob_data_over_parameter;
	mig1d_likelihood_list(i_propose) = mig1d_new_prob_data_over_parameter;
	sum_likelihood_list(i_propose) = sum_likelihood_new;
	prior_prob_list(i_propose) = prior_prob_new;
	posterior_prob_list(i_propose) = posterior_prob_new;
	running_t_list(i_propose) = toc;
	parameter_updated_list{i_propose} = parameter_updated_label;
	
	for fdn = fieldnames(param_list)'
		param_list.(fdn{1})(i_param_list) = param.(fdn{1});
	end
	i_param_list = i_param_list + 1;
	
	% screen output
	fprintf('%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%1.2e\t%1.2e\n', ...
		i_propose, ...
		accept, ...
		-sum_likelihood * error_tol^2, ...
		sum_likelihood, ...
		prior_prob, ...
		posterior_prob, ...
		parameter_updated_label, ...
		old_val, new_val);
	toc
	
	if (i_propose < 100 && mod(i_propose, n_propose_write) == 0) ...
			|| (i_propose >= 100 && mod(i_propose, 100*n_propose_write) == 0) ...
			|| i_propose == n_propose
		
		save(outfilepath ...
			, 'param_list' ...
			, 'wt_likelihood_list' ...
			, 'gal80d_likelihood_list' ...
			, 'mig1d_likelihood_list' ...
			, 'sum_likelihood_list'	...
			, 'prior_prob_list' ...
			, 'posterior_prob_list' ...
			, 'accept_list' ...
			, 'param_map' ...
			, 'param_prob_map' ...
			, 'i_propose' ...
			, 'running_t_list' ...
            , 'parameter_updated_list' ...
			, '-append' ...
			);
	end
end

end

function res = get_param_prob(param, parameter_update)
res = 0;
for i = 1:height(parameter_update)
	res = res + log( ...
			pdf('lognormal', param.(parameter_update{i, 'parameter_name'}{1}), ...
			log(parameter_update{i, 'prior_mean'}), ...
			parameter_update{i, 'prior_sigma'}) ...
            );
end
end

function prob = get_prob_data_over_parameter(param, trait, error_tol, fit_type)

output = evalGalPathway(param, trait, fit_type);
prob = -output.G1obj / error_tol^2;

end