function badsGAL(folder_name, jobtag, array_id, fit_type)

%% add path with subfolders to the search path
addpath(genpath('bads-master'))

%% load configuration file and GAL3, GAL4 level data
filepath = fullfile(folder_name, [jobtag, sprintf('_%03d',str2num(array_id)), '.mat']);
load(filepath)
GAL1_trait = trait;

load('../metaData/trait_extraction/GAL3pr_all_data.mat')
GAL3level = trait.basal_level;
load('../metaData/trait_extraction/GAL4pr_all_data.mat')
GAL4level = trait.basal_level;

%% setup parameters that to optimize, and custom output function
n_parameter = height(parameter_update);
parameter_name = parameter_update.parameter_name;
% get the starting point for those parameters
for i=1:n_parameter
    x0(i) = param_init.(parameter_name{i}); % be careful, DON'T use prior_mean as x0, because they are all the same values
end
prior_mean = parameter_update.prior_mean;
prior_sigma = parameter_update.prior_sigma;

% lb = zeros(1,n_parameter);  % the lower hard bounds are zeros
ub = 1000 .* prior_mean';
% ub(end-4:end) = 10;         % hill coefficient
% use 3 sigma interval as plausible boundary for each parameter
plb = prior_mean' ./ (exp(3 .* prior_sigma))';   % or exp(log(prior_mean) - 3 .* prior_sigma)'
pub = prior_mean' .* (exp(3 .* prior_sigma))';   % or exp(log(prior_mean) + 3 .* prior_sigma)'

% custom OPTIONS
opt = bads('defaults');     % get a default OPTIONS struct
opt.MaxFunEvals = n_propose;

% set outfile path
% if ~isdir('../results/badsOptim/RenanData/')
%     mkdir('../results/badsOptim/RenanData/');
% end
outfilepath = fullfile(...
    sprintf(...
    '%s-%s-%s.txt', ...
    jobtag, sprintf('%03d',str2num(array_id)), ...
    datestr(now, 'yymmdd_hh:MM') ...
    ) ...
    );

% custom output function
opt.OutputFcn = @(x, optimState, state) bads_outfun(x, optimState, state, outfilepath);

%% call bads to carry out optimization
bads(...
    @(x) mytarget(...
    update_param( param_init, parameter_name, x )...
    , GAL1_trait, 0, GAL3level, GAL4level, fit_type ) ...
    , x0, lb, ub, plb, pub ...
    , opt)

end

function sum_obj = mytarget(param, GAL1_trait, GAL2level, GAL3level, GAL4level, fit_type)

output = evalGalPathway_GAL34_changedR(param, GAL1_trait, GAL2level, GAL3level, GAL4level, fit_type);
sum_obj = output.sum_obj;

end
