function param = update_param( base_param, fields, values )

n_param_to_update = length(fields);
param = base_param;
for i_param_to_update = 1:n_param_to_update
    param.(fields{i_param_to_update}) = values(i_param_to_update);
end

end