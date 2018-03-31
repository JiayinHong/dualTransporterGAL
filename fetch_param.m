function [ param ] = fetch_param( param_base, param_id, param_list )
%   this function is called by 'supp_fig1.m', it takes in a base set of
%   paramters, and update the values taking from the param_list, indexed by
%   param_id
param_names = fieldnames(param_list);
for i = 1:length(param_names)
    param_field = param_names{i};
    param_base.(param_field) = param_list.(param_field)(param_id);
end
param = param_base;

end

