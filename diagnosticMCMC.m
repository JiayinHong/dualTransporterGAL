% this script is used to diagnostic how MCMC works
% 2017.07.13 by JH

%% fit GAL1, GAL3 & GAL4, long chains - 1000,000 iterations
mcmc_data_folder = '../results/GAL234-longChains/';
jobtags = {'wildtype_1c', 'wildtype_1r', 'wildtype_1r1c'...
            , 'MAP-wildtype_1c', 'MAP-wildtype_1r', 'MAP-wildtype_1r1c'...
            };
mcmc_result = load_mcmc_result(mcmc_data_folder, jobtags);

%% define a threshold for good fitting
load(mcmc_result{1,'filepath'}{1}, 'error_tol');
thresh = - 0.1 / error_tol ^2;  % i.e. obj < 0.1
filtered_result = mcmc_result(mcmc_result.param_prob_map > thresh, :);

n_example = height(filtered_result);
n_row = floor(n_example ^0.5);
n_col = ceil(n_example / n_row);

%% plot obj variation
figure
set(gcf, 'position', [424 136 1316 842])
for i_example = 1:n_example
    max_iter = filtered_result{i_example, 'max_iter'};
%     accept_list = filtered_result{i_example, 'accept_list'}{1};
%     average_accept_list = nan(max_iter,1);
%     for i_iter = 1:max_iter
%         average_accept_list(i_iter) = mean(accept_list(1:i_iter));
%     end
    param_prob_list = filtered_result{i_example, 'param_prob_list'}{1}(1:max_iter);
    obj_list = - param_prob_list .* error_tol ^2;

    subplot(n_row, n_col, i_example)
%     yyaxis left
    plot(1:100:max_iter, obj_list(1:100:max_iter), 'b.')
    xlabel('iteration')
    ylabel('obj')
    title(changeunderscore(filtered_result{i_example, 'jobtag'}{1}))
    xlim([0 max_iter])
    set(gca, 'FontSize', 12)
%     yyaxis right
%     plot(1:max_iter, average_accept_list)
%     ylabel('average_accept')

end


%% plot parameter value variation
figure
set(gcf, 'position', [462 83 1204 872])

i_example = 1;

max_iter = filtered_result{i_example, 'max_iter'};
load(filtered_result{i_example, 'filepath'}{1})
n_param = height(parameter_update);
update_parameters = parameter_update{:, 'parameter_name'};
n_row = floor(n_param ^0.5);
n_col = ceil(n_param / n_row);

for i_param = 1:n_param
    
    param_name = update_parameters{i_param};
    param_value = param_list.(param_name);
    param_value = param_value(1:max_iter);
    subplot(n_row, n_col, i_param)
    plot(1:max_iter, param_value, '.')
    xlabel('iteration')
    ylabel('param value')
    title(param_name)
    xlim([0 max_iter])
    set(gca, 'FontSize', 12)
    
end










