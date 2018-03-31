% this script is used to carry out single parameter scan and joint
% parameter scan for gal80d & mig1d mutants, starting from a 'best-fit'
% parameter value of wildtype, trying to figure out if it's possible to fit
% mutants phenotype by solely tuning one or two parameter values

% this script is of limited reusable value, consider to delete

% load trait tables
load('../metaData/trait_extraction/S288C-double_gradient/wildtype_all_data.mat')
wt_trait = trait;
load('../metaData/trait_extraction/S288C-double_gradient/gal80d_all_data.mat')
gal80d_trait = trait;
load('../metaData/trait_extraction/S288C-double_gradient/mig1d_all_data.mat')
mig1d_trait = trait;

clear trait
load(mcmc_result{1, 'filepath'}{1}, 'parameter_update');
param_names = parameter_update.parameter_name;

%% single parameter value scan
base_param = mcmc_result{17, 'param_map'};  % the best fit for wildtype

% set the output directory
mig1dpath = '../paramScanForMig1d/';
if ~isdir(mig1dpath)
    mkdir(mig1dpath)
end

gal80dpath = '../paramScanForGal80d/';
if ~isdir(gal80dpath)
    mkdir(gal80dpath)
end

% single scan
for param_name = param_names'
    ParamScan(base_param,param_name,mig1d_trait,mig1dpath);
    ParamScan(base_param,param_name,gal80d_trait,gal80dpath);
end

%% joint parameter value scan
base_param = mcmc_result{17, 'param_map'};

jointMig1dPath = '../jointScanForMig1d/';
if ~isdir(jointMig1dPath)
    mkdir(jointMig1dPath)
end

jointGal80dPath = '../jointScanForGal80d/';
if ~isdir(jointGal80dPath)
    mkdir(jointGal80dPath)
end

% hard-coded which joint parameters to scan
ParamScan(base_param, {'a80','ag80'}, gal80d_trait,jointGal80dPath)
ParamScan(base_param, {'a80','kf84'}, gal80d_trait,jointGal80dPath)
ParamScan(base_param, {'ag80','kf84'}, gal80d_trait,jointGal80dPath)
ParamScan(base_param, {'a1','ag1'}, gal80d_trait,jointGal80dPath)
ParamScan(base_param, {'a3','ag3'}, gal80d_trait,jointGal80dPath)

ParamScan(base_param, {'aR','KR1'}, mig1d_trait,jointMig1dPath)
ParamScan(base_param, {'aR','KRs'}, mig1d_trait,jointMig1dPath)
ParamScan(base_param, {'KR1','KRs'}, mig1d_trait,jointMig1dPath)

