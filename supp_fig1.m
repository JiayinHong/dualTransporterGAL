function MaxAndIndex = supp_fig1( mcmc_result, jobtag, trancation )
	% this function is used to plot 1) how the probability density varied along with the elongation of the chain;
	% 2) the histogram showing the iteration that MAP lies in; 3) pick up several chains to show how the 'MAP'
	% within 5k iterations, 10k iterations, 20k iterations and of the whole chain fitting the experiment data;
	% 4) violin plot shows the distribution of each parameter of MAP, of all the chains.

	%	data preprocess
	sorted_mcmc_result = sortrows(mcmc_result, 'param_prob_map', 'descend');
	mcmc_tmp = sorted_mcmc_result(ismember(sorted_mcmc_result.jobtag, jobtag), :);
	tmp_long_chain = mcmc_tmp(mcmc_tmp.max_iter >= trancation, :);
	param_prob_list_all = tmp_long_chain{1:end, 'param_prob_list'};
	myarray = [];

	n_chains = height(tmp_long_chain);
	
	Mall_list = nan(n_chains, 1);
	Iall_list = nan(n_chains, 1);
	M5k_list = nan(n_chains, 1);
	I5k_list = nan(n_chains, 1);
	M10k_list = nan(n_chains, 1);
	I10k_list = nan(n_chains, 1);
	M20k_list = nan(n_chains, 1);
	I20k_list = nan(n_chains, 1);

	for i_chain = 1:n_chains
		myarray(i_chain, :) = param_prob_list_all{i_chain};

		[Mall, Iall] = max(myarray(i_chain, :));	% the MAP of the whole chain
		Mall_list(i_chain, 1) = Mall;
		Iall_list(i_chain, 1) = Iall;

		[M5k, I5k] = max(myarray(i_chain, 1:5000));		% the MAP within the first 5000 iterations
		M5k_list(i_chain, 1) = M5k;
		I5k_list(i_chain, 1) = I5k;

		[M10k, I10k] = max(myarray(i_chain, 1:10000));
		M10k_list(i_chain, 1) = M10k;
		I10k_list(i_chain, 1) = I10k;

		[M20k, I20k] = max(myarray(i_chain, 1:20000));
		M20k_list(i_chain, 1) = M20k;
		I20k_list(i_chain, 1) = I20k;

	end

	MaxAndIndex = table(Mall_list, Iall_list, M5k_list, I5k_list, M10k_list, I10k_list, M20k_list, I20k_list, 'VariableNames', {'MAP_all', 'ind_all', 'MAP_5k', 'ind_5k', 'MAP_10k', 'ind_10k', 'Map_20k', 'ind_20k'});

	CT = cbrewer('qual', 'Set3', 8);	% generate a color palette
	
	%	shaded error bar plot of log probability density
	figure
	shadedErrorBar(1:size(myarray, 2), mean(myarray, 1), std(myarray), {'color', CT(1,:)});
	ylim([-50 0]);
	xlabel('iteration');
	ylabel('log\_probability\_density');
	title(jobtag);

	%	histogram showing the iteration that MAP lies in
	figure
	hist(Iall_list, 20)		% set 20 bins as default
	h = findobj(gca, 'Type', 'patch');
	h.Facecolor = CT(6,:);
	h.Edgecolor = 'w';
	xlabel('iteration that MAP lies in');
	ylabel('counts');
	title(jobtag);

	%	fetch parameter
	n_example = 4;	% show how many chains as examples

	param_map_list = cell(n_chains, 1);		% store all the parameters of map of all the chains into a cell array
	param_5k_list = cell(n_chains, 1);		% store all the parameters of best fitting within 5k iterations
	param_10k_list = cell(n_chains, 1);
	param_20k_list = cell(n_chains, 1);

	for i_chain = 1:n_chains
		filepath = tmp_long_chain{i_chain, 'filepath'}{1};
		load(filepath)
		ind_5k = I5k_list(i_chain);
		ind_10k = I10k_list(i_chain);
		ind_20k = I20k_list(i_chain);
		param_map = tmp_long_chain{i_chain, 'param_map'};
		param_5k = fetch_param(param_map, ind_5k, param_list);
		param_10k = fetch_param(param_map, ind_10k, param_list);
		param_20k = fetch_param(param_map, ind_20k, param_list);

		param_map_list{i_chain} = param_map;
		param_5k_list{i_chain} = param_5k;
		param_10k_list{i_chain} = param_10k;
		param_20k_list{i_chain} = param_20k;

	end

	%	MAP fitting
	figure
	set(gcf, 'position', [680   565   840   413])
	for i_example = 1:n_example

		subplot( n_example, 4, 4*(i_example-1)+1 )
		param_fitting_plot(param_5k_list{i_example}, trait);
		subplot( n_example, 4, 4*(i_example-1)+2 )
		param_fitting_plot(param_10k_list{i_example}, trait);
		subplot( n_example, 4, 4*(i_example-1)+3 )
		param_fitting_plot(param_20k_list{i_example}, trait);
		subplot( n_example, 4, 4*(i_example-1)+4 )
		param_fitting_plot(param_map_list{i_example}, trait);

	end

	%	violin plot, the parameter value is in log10 scale
	figure
	parameter_update = readtable('MCMC_parameter_config_including_all_n_d.csv');
	update_parameters = parameter_update{:, 'parameter_name'};
	param_map = tmp_long_chain{1, 'param_map'};	% just randomly pick a set of parameters as initial value
	all_param_value_map = param_map;	% make a struct to store all the parameter values in the 
										% corresponding field, the value will be overwritten later

	for i_param = 1:length(update_parameters)
		field_name = update_parameters{i_param};
		value_list = [];	% a value list to store all the values of a single parameter across all the chains
		for i_chain = 1:n_chains
			value_i_chain = log10(param_map_list{i_chain}.(field_name));
			value_list(i_chain) = value_i_chain;
		end
		all_param_value_map.(field_name) = value_list;
	end

	[~, catnames] = violinplot(all_param_value_map, 'ViolinColor', CT(7,:));
	%	this will make all the violin in the same color, default is a nice color cycle
	catnames = catnames(:);	% this is very important, be careful! the order of updated parameters
	% are not the same as the order of parameter fieldnames, so we need to store the category names

	%	plot the parameter value distribution of MAP within 10k
	hold on

	for i_param = 1:length(catnames)
		field_name = catnames{i_param};
		value_list = [];
		for i_chain = 1:n_chains
			value_i_chain = log10(param_10k_list{i_chain}.(field_name));
			value_list(i_chain) = value_i_chain;
		end

		plot(i_param, value_list, 'kd');
		hold on
	end

end