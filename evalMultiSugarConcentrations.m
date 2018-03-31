function y_ss_list = evalMultiSugarConcentrations( param, init_val, gluc_condition, galc_condition )

n_var = length(init_val);
n_condition = length(gluc_condition);
y_ss_list = nan(n_condition, n_var);

load_global;
opt = odeset('NonNegative',1:12);
accurate_thresh = 10^-8;
% opt = [];

for i_condition = 1:n_condition
    
    param.exglu = gluc_condition(i_condition) * perc_to_nm;
    param.exgal = galc_condition(i_condition) * perc_to_nm;
    odefunc = @(t,y)GALode5(t,y,param);
    curInitVal = init_val;
    
    if param.exglu == 0
        % R*=0, glu=0
%         tmp = [ones(1,7),0,1,1,0,1];
        tmp = [ones(1,9),0,1,1];
        curInitVal = curInitVal .* tmp;
    end
    if param.exgal == 0
        % Gal3*=0, C83=0, gal=0
%         tmp = [ones(1,5),0,1,1,0,1,1,0];
        tmp = [ones(1,5),0,1,0,1,1,0,1];
        curInitVal = curInitVal .* tmp;
    end
    if param.aR == 0
        % Repressor=0, R*=0
%         tmp = [ones(1,6),0,0,ones(1,4)];
        tmp = [ones(1,6),0,ones(1,5)];
        curInitVal = curInitVal .* tmp;
    end
    if param.a80 == 0 && param.ag80 == 0
        % Gal80=0, C83=0, C84=0
%         tmp = [ones(1,4),0,ones(1,3),0,0,1,1];
        tmp = [ones(1,4),0,ones(1,2),0,0,1,1,1];
        curInitVal = curInitVal .* tmp;
    end
    
    [~, y_current] = ode15s(odefunc, [0 10000], curInitVal, opt);
    y_current(y_current<accurate_thresh) = 0;   % omit values that are too small
    y_ss_list(i_condition,:) = y_current(end,:);
    
    %     % the previous version of ode solver, consider to delete
    %     option = 'expt_setup';
    %     y_current = run_till_steady_state( @(t, y) GAL_sec(t, y, param), curInitVal, option );
    %     y_ss_list(i_condition, :) = y_current;
end

end

