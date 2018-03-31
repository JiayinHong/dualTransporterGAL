function [ y_ss ] = run_till_steady_state( odefunc, initial_value, option )

time_span_increase = 1000;  % further increase if not at ss.
i_max_iteration = 0;        % current iteration
n_max_iteration = 1;        % maximum iterations

if nargin < 3
    fprintf('run for 8 hours!\n')
    time_span_init = 60 * 8;  % initial try
    
else
    switch option
        case 'initial_run'
            time_span_init = 60 * 80;
        case 'expt_setup'
            time_span_init = 60 * 8;
    end
end



done = false;
t_current = 0;
y_current = initial_value;
t_target = time_span_init;

while ~done && i_max_iteration < n_max_iteration
    i_max_iteration = i_max_iteration+1;
    
    opt = odeset('NonNegative', 1:12);
    [t, y] = ode15s( odefunc, [t_current, t_target], y_current, opt);
    
    done = check_stable(y(end-1, :), y(end, :));
    
    t_current = t(end);
    y_current = y(end, :);
    
    t_target = t_current + time_span_increase;
end

y_ss = y_current;
% if ~done
%     warning('not converging')
%     y_ss = y_current;
%     return
% else
%     y_ss = y_current;
% end

end