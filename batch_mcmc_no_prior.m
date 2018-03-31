function batch_mcmc_no_prior( folder_name, jobtag, array_id, fit_type )

task_id = str2double(array_id);
filepath = fullfile(folder_name, [jobtag, num2str(task_id, '_%03d'), '.mat']);
load(filepath);

mcmc_for_GAL234( trait, param_init, parameter_update, fit_type, 'n_propose', n_propose, 'jobtag', jobtag, 'arrayid', array_id );

end


