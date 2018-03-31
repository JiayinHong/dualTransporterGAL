function batch_fminsearch_on_slurm( folder_name, jobtag, array_id, fit_type )

    task_id = str2double(array_id);
    filepath = fullfile(folder_name, [jobtag, num2str(task_id, '_%03d'), '.mat']);
    load(filepath);
    fminsearchGAL( trait, param_init, parameter_update, jobtag, array_id, fit_type );
    
end
