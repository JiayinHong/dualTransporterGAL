function singleParamScan(arrayID)
% this function is called directly from the shell script of the same name,
% for a batch of strains to implement single parameter scan simultaneously
addpath(genpath('../../matlab_script/'))

strainID = str2double(arrayID);
clear param
load('best_param.mat','param');
param_names = fieldnames(param);
singleParamScanHelper(param,param_names,strainID)

end

function singleParamScanHelper(base_param,param_names,strainID)
% calculate single parameter scan obj values and store
traitsFolder = '../../selectedTraits/';
allTraits = dir(traitsFolder');
allTraits = allTraits([allTraits.isdir]~=1);
obj = struct();     % store the results in a struct, with a fieldname naming after the strainID
clear trait
fileName = allTraits(strainID).name;
load(fullfile(traitsFolder, fileName),'trait')     % load certain strain trait

% logiArray = cellfun(@(x) regexp(x,'n*'), param_names, 'UniformOutput', false);
% ind = find(cellfun(@isempty,logiArray));    % to find parameters that are not hill coefficients

% for i_param = ind'
for i_param = 1:length(param_names)
    tmp(i_param,:) = ParamScan2(base_param, param_names(i_param), trait);
end
fieldname = sprintf('str%02d',strainID);
obj.(fieldname) = tmp;

% make the directory and save results
% foldername = '../paramScanResults/';
% if ~isdir(foldername)
%     mkdir(foldername)
% end
% save(fullfile(foldername, sprintf('res-%02d.mat',strainID)),'obj','allTraits')
save(sprintf('res-%02d.mat',strainID),'obj','allTraits')

% the reason to save the variable 'allTraits' is that the sequence of the
% strain names here is very important, since I am not using a string to
% name each strain, which would cause invalid naming. Here all the names
% are like 'str1', 'str2', etc. Hence, we need to save 'allTraits' to make
% sure the sequence of these strains is not messy

end


