function shell_on_slurm_generator( algorithm, jobtag, folder_name, fit_type )
%   this function is called by 'generate_shell_on_slurm'

str1 = '#!/bin/bash\n';
if algorithm == 'bads'
    str2 = '#SBATCH -p short\n';
    str3 = '#SBATCH -t 0-12:00\n';
else
    str2 = '#SBATCH -p medium\n';
    str3 = '#SBATCH -t 5-00:00\n';
end
str4 = '#SBATCH -c 1\n';
str5 = '#SBATCH --mem-per-cpu=10G\n';
str6 = '#SBATCH -e ../err_files/%%A_%%a.err\n';
str7 = '#SBATCH -o ../out_files/%%A_%%a.out\n';
str8 = 'module load gcc/6.2.0 matlab/2017b\n';
long_string = [str1, str2, str3, str4, str5, str6, str7, str8];

% to avoid the conflict of shell scripts name
script_name = sprintf('slurm_%s_%s.sh', algorithm, jobtag);
% script_name = sprintf('slurm_%s_%s_withPrior.sh', algorithm, jobtag);             % case 'with Prior'
fid = fopen(fullfile('../shell_script',script_name), 'w');
fprintf(fid, long_string);
prefix = 'matlab -nodesktop -nosplash -nojit -r ';
switch algorithm
    case 'mcmc'
%         fprintf(fid, '%s "batch_mcmc_no_prior ''%s'' ''%s'' ''${SLURM_ARRAY_TASK_ID}'' ''%s'' "', prefix, folder_name, jobtag, fit_type);

        fprintf(fid, '%s "batch_mcmc_on_slurm ''%s'' ''%s'' ''${SLURM_ARRAY_TASK_ID}'' ''%s'' "', prefix, folder_name, jobtag, fit_type);
    case 'bads'
        fprintf(fid, '%s "badsGAL ''%s'' ''%s'' ''${SLURM_ARRAY_TASK_ID}'' ''%s'' " ', prefix, folder_name, jobtag, fit_type);
    case 'fminsearch'
        fprintf(fid, '%s "batch_fminsearch_on_slurm ''%s'' ''%s'' ''${SLURM_ARRAY_TASK_ID}'' ''%s'' "', prefix, folder_name, jobtag, fit_type);
end

end

