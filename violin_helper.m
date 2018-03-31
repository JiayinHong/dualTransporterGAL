function [ log10_param_value_list, linear_param_value_list ] = violin_helper( param_map_list )
% this function is to help making separate violin plot by scale, it defines what parameters you
% want to plot in log scale and what in linear scale, and generate two corresponding parameter
% value list that can served as the input to function 'violinplot'

log10_scale_param = {'a1', 'a3', 'a4', 'a80', 'aR', 'ag1', 'ag3', 'ag4', ...
    'ag80', 'KMglu', 'KMgal', 'kf3', 'kf83', 'kf84', 'kfR', 'KG1', 'KG3', ...
    'KG80', 'KR1', 'KR3', 'KR4'};
linear_scale_param = {'d', 'n1', 'n2', 'n3', 'n80', 'nR1', 'nR3', 'nR4'};

log10_param_value_list = struct();
linear_param_value_list = struct();

n_chains = length(param_map_list);
for field_name = log10_scale_param
    value_list = [];	% a value list to store all the values of a single parameter across all the chains
    for i_chain = 1:n_chains
        value_i_chain = log10(param_map_list{i_chain}.(field_name{1}));
        value_list(i_chain) = value_i_chain;
    end
    log10_param_value_list.(field_name{1}) = value_list;
end

for field_name = linear_scale_param
    value_list = [];
    for i_chain = 1:n_chains
        value_i_chain = param_map_list{i_chain}.(field_name{1});
        value_list(i_chain) = value_i_chain;
    end
    linear_param_value_list.(field_name{1}) = value_list;
    
end

