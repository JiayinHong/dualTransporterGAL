function mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags)
% this version is only suitable for the special version of mcmc that
% without prior

% load all result files

mcmc_result = table();
files = dir( fullfile(mcmc_data_folder, '*.mat') );

% filename_pat = '([\d_\w]+)_\d{8}_\d{6}_\d{3}_\d+.mat';
% example: 'MCMC_BC187_random_init_0217_20170217_161008_731_62002.mat'

% filename_pat = '([\d_\w]+)-\d{1,3}_\d{8}_\d{6}.mat';
% example: 'MCMC_multiple_trait-11_20170410_155138.mat'

% filename_pat = '([\d_\w]+)-\d{1,3}_\d{8}_\d{2}:\d{2}.mat';
% example: 'MCMC_wt_1c-13_20170422_20:53.mat'

% filename_pat = '([\d_\w]+)-\d{3}-\d{6}_\d{2}:\d{2}.mat';
% example: 'wildtype_1r1c-001-170502_18:12.mat'

filename_pat = '([\w\d_\w]+)-\d{3}-\d{6}_\d{2}.*\d{2}.mat';
% example: 'sequestrate-wildtype_1r-024-170802_19:57.mat'


for i_file = 1:length(files)
    filename = files(i_file).name;
    filepath = fullfile(mcmc_data_folder, filename);
    tok = regexp(filename, filename_pat, 'tokens');
    if isempty(tok)
        fprintf('File name format warning: %s\n', filename)
    else
        if ismember(tok{1}, jobtags)
            try
                load(filepath);
            catch
                warning('Not able to load %s\n', filepath);
                continue
            end
            
            max_iter = i_propose;
            average_accept = nanmean(accept_list);
            
            if strcmp(mcmc_version, '3.0')
%                 [max_likelihood,IndMAP] = max(sum_likelihood_list);
                mcmc_result = [ ...
                    mcmc_result; ...
                    table( ...
                    {filepath}, max_iter, param_map, param_prob_map, {jobtag}, average_accept, {trait}, ...
                    'VariableNames', {'filepath', 'max_iter', 'param_map', 'param_prob_map', 'jobtag', 'average_accept', 'trait'} ...
                    ) ...
                    ];
            
            elseif regexp(mcmc_version, '\w*prior_included$')
                [~,index] = max(param_prob_list);
                map_data_over_param = prob_data_over_parameter_list(index);
                map_param_over_prior = prob_parameter_over_prior_list(index);
                mcmc_result = [ ...
                    mcmc_result; ...
                    table( ...
                    {filepath}, max_iter, param_map, param_prob_map, map_data_over_param, map_param_over_prior, {jobtag}, {prob_data_over_parameter_list}, average_accept, param_list, ...
                    'VariableNames', {'filepath', 'max_iter', 'param_map', 'param_prob_map', 'map_data_over_param', 'map_param_over_prior', 'jobtag', 'prob_data_over_parameter_list', 'average_accept', 'param_list'} ...
                    )...
                    ];
                
%             elseif strcmp(mcmc_version, 'Simulated_Annealing')
            elseif strcmp(mcmc_version, 'AlsoFitGAL34')
                mcmc_result = [ ...
                    mcmc_result; ...
                    table( ...
                    {filepath}, max_iter, param_map, param_prob_map, {jobtag}, {prob_data_over_parameter_list}, {accept_list}, average_accept, {GAL1_trait}, ...
                    'VariableNames', {'filepath', 'max_iter', 'param_map', 'param_prob_map', 'jobtag', 'param_prob_list', 'accept_list', 'average_accept', 'trait'} ...
                    )...
                    ];
            
            else
                [~, i] = max(prob_data_over_parameter_list);
                try
                    simulation_result_linear_map = simulation_result_linear_list(i);
                catch
                    simulation_result_linear_map = nan;
                end
                
                mcmc_result = [...
                    mcmc_result; ...
                    table(...
                    {filepath}, max_iter, param_map, param_prob_map, {jobtag}, {prob_data_over_parameter_list}, average_accept, {trait}, simulation_result_linear_map, ...
                    'VariableNames', {'filepath', 'max_iter', 'param_map', 'param_prob_map', 'jobtag', 'param_prob_list', 'average_accept', 'trait', 'simulation_result_linear_map'}...
                    )...
                    ];
            end
        else
            % pass
        end
    end
end

fprintf('done!\n')
