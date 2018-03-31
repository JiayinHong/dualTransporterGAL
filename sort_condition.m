function [ gluc_list, galc_list, id ] = sort_condition( gluc_list, galc_list )
% this function is called by 'eval_param'
% this function sort the sugar condition as the ratio of galactose/glucose ascending

[~,id] = sort(galc_list./gluc_list, 'ascend');

galc_list = galc_list(id);
gluc_list = gluc_list(id);

end