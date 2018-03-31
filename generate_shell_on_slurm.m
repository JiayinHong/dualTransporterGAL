% this script is used to generate shell scripts to call the mcmc function
% to run on slurm
% the core of the script has been modified to a function, take folder name, jobtag and fit type as input

algorithm_list = {'mcmc', 'fminsearch'};
jobtag_list = {'wildtype_1r', 'mig1d_1r', 'gal80d_1r', 'triple_fit_1r', ...
    'wildtype_1c', 'mig1d_1c', 'gal80d_1c', 'triple_fit_1c', ...
    'wildtype_1r1c', 'mig1d_1r1c', 'gal80d_1r1c', 'triple_fit_1r1c'};
for algorithm = algorithm_list
    algorithm = algorithm{1};
    
    for jobtag = jobtag_list
        jobtag = jobtag{1};
        
        if regexp(jobtag, 'triple_fit')
            folder_name = '../metaData/mutant_and_wt_triple_fit';
        else
            folder_name = '../metaData/random_init_mutant_and_wt';
        end
        
        if regexp(jobtag, '\w*_1r$') % match words ending with 1r
            fit_type = 'one_row';
        elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
            fit_type = 'one_column';
        elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
            fit_type = 'one_cross';
        elseif regexp(jobtag, '\w*_96well$')    % match words ending with 96well
            fit_type = '96well';
        end
        
        shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
    end
end

%% generate shell scripts to fit GAL1,3,4 in all 96 well data
algorithm = 'mcmc';
% jobtag_list = {'wildtype_96well', 'mig1d_96well', 'gal80d_96well'};
jobtag_list = {'BC187_Kayla'};
folder_name = '../metaData/biTrans_addHXT_BC187/';
fit_type = '96well';
for jobtag = jobtag_list
    jobtag = jobtag{1};
shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
end

%% re-run MCMC to fit BC&YJM glucose gradient data
algorithm = 'mcmc';
jobtag_list = {'BC187', 'YJM978'};
folder_name = '../metaData/BCandYJM_single_grad';

for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    fit_type = 'gluc_gradient';
    
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
end

%% generate shell scripts for GAL1,3,4 fitting
algorithm = 'mcmc';
jobtag_list_varyN = {'varyN-wildtype_1r', 'varyN-wildtype_1c', 'varyN-wildtype_1r1c'};
jobtag_list_sequestrate = {'sequestrate-wildtype_1r', 'sequestrate-wildtype_1c', 'sequestrate-wildtype_1r1c'};

for jobtag = jobtag_list_varyN
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
    end
    
    folder_name = '../metaData/Aug1st-fitGAL134-vary-n/';
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
    
end

for jobtag = jobtag_list_sequestrate
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
    end
    
    folder_name = '../metaData/Aug1st-fitGAL134-pure-sequestration/';
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
end


%% use alpha*KMglu to replace KMgal, GAL1,3,4 fitting, vary hill coefficients
algorithm = 'mcmc';
% test different step size and the span of prior distribution
jobtag_list = {'small-wildtype_1r', 'small-wildtype_1c', 'small-wildtype_1r1c'...
              ,'medium-wildtype_1r', 'medium-wildtype_1c', 'medium-wildtype_1r1c'...
              ,'large-wildtype_1r', 'large-wildtype_1c', 'large-wildtype_1r1c'};

for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
    end
    
    folder_name = '../metaData/fitGAL134-TestStepSize/';
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
    
end

%% change the formula of Mig1
algorithm = 'mcmc';
jobtag_list = {'medium-wildtype_1r', 'medium-wildtype_1c', 'medium-wildtype_1r1c'...
               , 'medium-gal80d_1r', 'medium-gal80d_1c' ...
               , 'medium-mig1d_1r', 'medium-mig1d_1c'};
for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
    end
    
    folder_name = '../metaData/fitGAL134-changeRform/';
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
    
end

%% Use internal glucose to replace Mig1 in all the equations
algorithm = 'mcmc';
jobtag_list = {'medium-wildtype_1r', 'medium-wildtype_1c', 'medium-wildtype_1r1c'};
               
for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$') % match words ending with 1r
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$') % match words ending with 1c
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')   % match words ending with 1r1c
        fit_type = 'one_cross';
    end
    
    folder_name = '../metaData/biTrans_addHXT/';
    shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
    
end

%% test what if we don't force aR=0 or a80=0, ag80=0, using mcmc
jobtag_list = {'all_update_mig1d_1r', 'all_update_gal80d_1r', ...
    'all_update_mig1d_1c', 'all_update_gal80d_1c', ...
    'all_update_mig1d_1r1c', 'all_update_gal80d_1r1c'};
folder_name = '../metaData/random_init_mutant_and_wt';
for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$')
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$')
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')
        fit_type = 'one_cross';
    end
    
    shell_on_slurm_generator( 'mcmc', jobtag, folder_name, fit_type )
end

%% test what if we don't force aR=0 or a80=0, ag80=0, using fminsearch
jobtag_list = {'all_update_mig1d_1r', 'all_update_gal80d_1r', ...
    'all_update_mig1d_1c', 'all_update_gal80d_1c', ...
    'all_update_mig1d_1r1c', 'all_update_gal80d_1r1c'};
folder_name = '../metaData/random_init_mutant_and_wt';
for jobtag = jobtag_list
    jobtag = jobtag{1};
    
    if regexp(jobtag, '\w*_1r$')
        fit_type = 'one_row';
    elseif regexp(jobtag, '\w*_1c$')
        fit_type = 'one_column';
    elseif regexp(jobtag, '\w*_1r1c$')
        fit_type = 'one_cross';
    end
    
    shell_on_slurm_generator( 'fminsearch', jobtag, folder_name, fit_type )
end