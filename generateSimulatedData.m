% this script is used to generate a trait table using simulated data, so
% that there's no experimental fluctuation or noise. The next step is to
% add noise based on the Coefficient of Variation or Standard Deviation of
% experimental data.
% Craeted by JH on 01.05.2018

%% load a good fit of BC187
filepath = mcmc_result{1,'filepath'}{1};
load(filepath, 'param_map', 'parameter_update', 'GAL1_trait')
param_names = parameter_update.parameter_name;
base_param = param_map;
mask_basal = GAL1_trait.mask_basal;
mask_induction = GAL1_trait.mask_induction;
basal_frac = GAL1_trait.basal_frac;
ind_frac = GAL1_trait.ind_frac;

%% setup the sugar titration gradient
gluc_gradient = [2 .^ [0:-1:-6], 0]';
galc_gradient = [0, 2 .^ [-8:1:2]];
gluc = gluc_gradient * ones(1,12);
galc = ones(8,1) * galc_gradient;
gluc = gluc(:);
galc = galc(:);

%% generate none noisy trait table
output = getInitParamsGalPathway(base_param);
y_ss_Glu = evalMultiSugarConcentrations(base_param, output.y0_Glu, gluc, galc);
y_ss_Gal = evalMultiSugarConcentrations(base_param, output.y0_Gal, gluc(end:-1:1), galc(end:-1:1));  % due to hysteresis?
y_ss_Gal = y_ss_Gal(end:-1:1,:);
basal_level = y_ss_Glu(:,1);
ind_level = y_ss_Gal(:,1);

basal_level = log(basal_level ./3000) + 7.34;   % funcional inverse of logyfp_to_nm
ind_level = log(ind_level ./3000) + 7.34;

trait = table(basal_level, ind_level, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc ...
    , 'VariableNames', {'basal_level', 'ind_level', 'basal_frac', 'ind_frac', 'mask_basal', 'mask_induction' ...
    , 'gluc', 'galc'});

save(fullfile('../metaData/trait_extraction/', 'competitiveNoneNoisy.mat'), 'trait')

%% generate add noise trait table
load('../metaData/trait_extraction/BC187_Kayla_Nov29.mat')
load_global
ind_off = [trait.mask_basal] ~= 0;
linearLo = logyfp_to_nm([trait{ind_off,'basal_level'}]);
linearHi = logyfp_to_nm([trait{[48:8:96],'ind_level'}]);
stdOFF = std(linearLo);
stdHiGal = std(linearHi);
meanOFF = mean(linearLo);
meanHiGal = mean(linearHi);
cvOFF = stdOFF/meanOFF;         % in Kayla's BC187 data, cvOFF = 0.3410
cvHiGal = stdHiGal/meanHiGal;   % in Kayla's BC187 data, cvHiGal = 0.0523
% I chose to use the smaller variance as the noise, i.e. cvHiGal = 0.0523

basal_level = y_ss_Glu(:,1);
ind_level = y_ss_Gal(:,1);
noisy_basal_level = basal_level .* (1+cvHiGal*randn(length(basal_level),1));
noisy_ind_level = ind_level .* (1+cvHiGal*randn(length(ind_level),1));
% if noisy basal level is higher than induced level, then swap them as
% below
tmp = nan(length(basal_level),1);
id1 = find(noisy_basal_level>noisy_ind_level);
tmp(id1) = noisy_ind_level(id1);
noisy_ind_level(id1) = noisy_basal_level(id1);
noisy_basal_level(id1) = tmp(id1);

% convert to fluorescent intensity in arbitrary unit, and save the trait
% table
noisy_basal_level = log(noisy_basal_level ./3000) + 7.34;   % funcional inverse of logyfp_to_nm
noisy_ind_level = log(noisy_ind_level ./3000) + 7.34;

trait = table(noisy_basal_level, noisy_ind_level, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc ...
    , 'VariableNames', {'basal_level', 'ind_level', 'basal_frac', 'ind_frac', 'mask_basal', 'mask_induction' ...
    , 'gluc', 'galc'});

save(fullfile('../metaData/trait_extraction/', 'competitiveAddNoise.mat'), 'trait')



