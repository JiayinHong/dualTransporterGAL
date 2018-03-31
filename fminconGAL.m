function fminconGAL( trait, param_init, parameter_update, jobtag, array_id, fit_type )

parameter_name = parameter_update{:, 'parameter_name'};
param_init_value = nan(1,length(parameter_name));
for i_param = 1:length(parameter_name)
    param_name = parameter_name{i_param};
    param_init_value(i_param) = param_init.(param_name);
end

if ~isdir('../results/fmincon');
    mkdir('../results/fmincon');
end

task_id = str2double(array_id);
outfilepath = fullfile(...
    '../results/fmincon/', ...
    sprintf(...
    '%s-%s-%s.txt', ...
    jobtag, num2str(task_id, '%03d'), ...
    datestr(now, 'yymmdd_hh:MM') ...
    ) ...
    );

% lb = zeros(1, length(parameter_name));
% ub = [];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% nonlcon = [];

lb = zeros(1, length(parameter_name));
ub = [];
A = -eye(length(parameter_name));
b = zeros(length(parameter_name),1);
Aeq = [];
beq = [];
nonlcon = [];

opts = optimoptions(@fmincon ...
    , 'Display', 'iter-detailed'...
    , 'OutputFcn', @(x, optimValues, state) fmincon_outfun(x, optimValues, state, outfilepath));

fmincon(...
    @(x) myfun(...
    update_param( param_init, parameter_name, x )...
    , trait, fit_type ) ...
    , param_init_value ...
    , A, b ...      % if there is a linear inequality constraint
    , Aeq, beq ...  % if there is a linear equality constraint
    , lb, ub ...    % if there are bound constraints
    , nonlcon ...
    , opts)

end

function sum_obj = myfun(param, trait, fit_type)
output = evalGalPathway(param, trait, fit_type);
sum_obj = output.sum_obj;
end