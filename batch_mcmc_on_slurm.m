function batch_mcmc_on_slurm( folder_name, jobtag, array_id, fit_type )
addpath(genpath('../../matlab_script/'))
task_id = str2double(array_id);
filepath = fullfile(folder_name, [jobtag, num2str(task_id, '_%03d'), '.mat']);
load(filepath);

if regexp(folder_name, 'triple_fit')
    simulateGALPathway( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
else
    mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
end


% elseif strcmp(folder_name, '../metaData/fitGAL134-changeRform/')
%     mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/fitGAL134-removeR/')
%     mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/fitGAL134-TestStepSize/')
%     mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/Aug1st-fitGAL134-vary-n/')
%     mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/Aug1st-fitGAL134-pure-sequestration/')
%     mcmc_for_GAL234_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/random_init_mutant_and_wt')
%     %     mcmc_without_prior( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     mcmc_for_GAL234( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% elseif strcmp(folder_name, '../metaData/BCandYJM_single_grad')
%     mcmc_prior_included( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );
%     
% else error('there''s something wrong, check folder name\n')
%     
% end