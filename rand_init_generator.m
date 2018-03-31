function rand_init_generator( n_init, trait, n_propose, base_param, parameter_update, jobtag, folder )
%   use this function to generate random initial state for later use like
%   mcmc or fminsearch, given trait, base_param, parameter_update, jobtag etc. 
%   this function is called by 'submit_to_orchestra'

error_tol = .15;

if nargin == 7
    folder_random_init = folder;
else
    folder_random_init = '../metaData/random_init_mutant_and_wt';
end
    
if ~exist(folder_random_init)
    mkdir(folder_random_init)
end

for i_init = 1:n_init
    parameter_val = nan(height(parameter_update), 1);
    for i = 1:height(parameter_update)
        % add the following line to make sure even if I forget to change
        % the csv config file, the prior mean will be updated, based on the
        % base parameter value
        
%         parameter_update{i, 'prior_mean'} = base_param.(parameter_update{i, 'parameter_name'}{1});
        
        parameter_val(i) = random('lognormal', ...
            log(parameter_update{i, 'prior_mean'}), ...
            parameter_update{i, 'prior_sigma'} ...
            );
    end
    param_init = update_param(base_param, parameter_update.parameter_name, parameter_val);
    
    save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat'])...
        , 'parameter_update', 'param_init' ...
        , 'error_tol', 'n_propose' ...
        , 'trait' ...
        , 'jobtag'...
        )
end

% the following part is to check when sigma is equal to zero, if the
% generated value is what it is supposed to be

% for i_init = 1:n_init
%     parameter_val = nan(height(parameter_update), 1);
%     for i = 1:height(parameter_update)
%         parameter_val(i) = random('lognormal', ...
%             log(parameter_update{i, 'prior_mean'}), ...
%             0 ...
%             );
%     end
%     param_init = update_param(base_param, parameter_update.parameter_name, parameter_val);
%     
%     save(fullfile(folder_random_init, [jobtag, num2str(i_init, '_%03d'), '.mat'])...
%         , 'parameter_update', 'param_init' ...
%         , 'error_tol', 'n_propose' ...
%         , 'trait' ...
%         , 'jobtag'...
%         )
% end

end

