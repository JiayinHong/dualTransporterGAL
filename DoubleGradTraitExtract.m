%% load data

clear
load('../RawData/mig1d_gal80d_wt_data.mat')

% reverse the glucose gradient, so that from top to bottom, glucose
% decreases, while from left to right, galactose increases
mig1d = cell(8,12);
for i_row = 8:-1:1
    for j_col = 1:12
        tmp = {alldata(3,i_row,j_col,2).yfp};
        mig1d{9-i_row,j_col} = tmp;
    end
end
% the dimensions of 'alldata': first dimension-strain type, 3 for mig1d,
% 4 for gal80d, and 5 for wildtype; second dimension-glucose gradient;
% third dimension-galactose gradient; fourth dimension-1 for reference,
% 2 for query
gal80d = cell(8,12);
for i_row = 8:-1:1
    for j_col = 1:12
        tmp = {alldata(4,i_row,j_col,2).yfp};
        gal80d{9-i_row,j_col} = tmp;
    end
end

wildtype = cell(8,12);
for i_row = 8:-1:1
    for j_col = 1:12
        tmp = {alldata(5,i_row,j_col,2).yfp};
        wildtype{9-i_row,j_col} = tmp;
    end
end

%% make the trait table for double gradient condition, mig1 delete, gal80 delete, and wildtype

plate_type_list = {'mig1d', 'gal80d', 'wildtype'};
for plate_name = plate_type_list;
    plate_name = plate_name{1};
    
    % prep trait
    basal_level = nan(8,12);
    ind_level = nan(8,12);
    modality = nan(8,12);
    basal_frac = nan(8,12);
    ind_frac = nan(8,12);
    mask_basal = 0.5 .*ones(8,12);
    mask_induction = 0.5 .*ones(8,12);
    
    gluc_gradient = [2 .^[0:-1:-6], 0]';
    galc_gradient = [0, 2 .^[-8:1:2]];
    gluc = gluc_gradient * ones(1,12);
    galc = ones(8,1) * galc_gradient;
    
    [hist_bin_edge, hist_bin_center] = GetHistBin(0.05, 9.95, 0.1);
    
    switch plate_name
        
        case 'mig1d'
            
            cutoff_mig1d = nan(8,12);
            ref_data = mig1d{1,1}{1};   % fully OFF
            log_ref_data = log(ref_data);   % natural logarithm
            counts_total_ref = size(ref_data,1);    % the total counts of ref data
            counts_list_ref = histcounts(log_ref_data,hist_bin_edge);
            % get a list of counts in each bin
            normalized_counts_ref = counts_list_ref ./ counts_total_ref;
            threshold = mean(log_ref_data);
            % ignore if there are more counts (normalized) of query than ref
            % data beneath the threshold
            
            for i_row = 1:8
                for j_col = 1:12
                    
                    query_data = mig1d{i_row,j_col}{1};
                    log_query_data = log(query_data);   % natural logarithm
                    counts_total_query = size(query_data,1);
                    counts_list_query = histcounts(log_query_data,hist_bin_edge);
                    normalized_counts_query = counts_list_query ./ counts_total_query;
                    
                    delta_freq = normalized_counts_query - normalized_counts_ref;
                    % the difference between query and fully OFF, regarded as the ON fraction of query data
                    weight_on_peak = delta_freq;    % the weight of ON peak in each bin
                    weight_on_peak(weight_on_peak<0) = 0;
                    weight_on_peak(weight_on_peak>0 & ...
                        [1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
                    % 0.1 is the bin size
                    i_cutoff = find(weight_on_peak>0,1);
                    % find in which bin is the first ON fraction lies in
                    if i_cutoff
                        cutoff_mig1d(i_row,j_col) = hist_bin_edge(i_cutoff);
                    end
                    % get the cutoff line position
                    ind_frac(i_row,j_col) = sum(weight_on_peak);
                    ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
                    basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
                    weight_off_peak = normalized_counts_query - weight_on_peak;
                    basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);
                    
                end
            end
            
            detect_limit = basal_frac(8,12);
            % the lower right corner well should be fully ON, where
            % basal_frac should be zero
            % the (8,10) well is hard-coded, since it should be fully ON
            for i_row = 1:8
                for j_col = 1:12
                    if basal_frac(i_row,j_col) <= max(basal_frac(8,10), detect_limit)
                        mask_basal(i_row,j_col) = 0;
                        mask_induction(i_row,j_col) = 1;
                    elseif ind_frac(i_row,j_col) <= max(basal_frac(8,10), detect_limit)
                        mask_induction(i_row,j_col) = 0;
                        mask_basal(i_row,j_col) = 1;
                    end
                end
            end
            
            
        case 'gal80d'
            
            cutoff_gal80d = nan(8,12);
            % use the upper left corner of the wildtype data as reference
            ref_data = wildtype{1,1}{1};
            log_ref_data = log(ref_data);
            counts_list_ref = histcounts(log_ref_data,hist_bin_edge);
            counts_total_ref = size(ref_data,1);
            normalized_counts_ref = counts_list_ref ./ counts_total_ref;
            threshold = mean(log_ref_data);
            
            for i_row = 1:8
                for j_col = 1:12
                    
                    query_data = gal80d{i_row,j_col}{1};
                    log_query_data = log(query_data);
                    counts_list_query = histcounts(log_query_data,hist_bin_edge);
                    counts_total_query = size(query_data,1);
                    normalized_counts_query = counts_list_query ./ counts_total_query;
                    delta_freq = normalized_counts_query - normalized_counts_ref;
                    weight_on_peak = delta_freq;
                    weight_on_peak(weight_on_peak<0) = 0;
                    weight_on_peak(weight_on_peak>0 & ...
                        [1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
                    i_cutoff = find(weight_on_peak>0,1);
                    
                    if i_cutoff
                        cutoff_gal80d(i_row,j_col) = hist_bin_edge(i_cutoff);
                    end
                    
                    ind_frac(i_row,j_col) = sum(weight_on_peak);
                    ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
                    basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
                    weight_off_peak = normalized_counts_query - weight_on_peak;
                    basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);
                    
                end
            end
            
            detect_limit = basal_frac(8,12);
            for i_row = 1:8
                for j_col = 1:12
                    if basal_frac(i_row,j_col) <= max(basal_frac(1,12), detect_limit)
                        mask_basal(i_row,j_col) = 0;
                        mask_induction(i_row,j_col) = 1;
                    elseif ind_frac(i_row,j_col) <= max(basal_frac(1,12), detect_limit)
                        mask_induction(i_row,j_col) = 0;
                        mask_basal(i_row,j_col) = 1;
                    end
                end
            end
            
            
        case 'wildtype'
            
            cutoff_wildtype = nan(8,12);
            ref_data = wildtype{1,1}{1};
            log_ref_data = log(ref_data);
            counts_list_ref = histcounts(log_ref_data,hist_bin_edge);
            counts_total_ref = size(ref_data,1);
            normalized_counts_ref = counts_list_ref ./ counts_total_ref;
            threshold = mean(log_ref_data);
            % if delta freqency in bins beneath threshold > 0,
            % drop it, make it equal to 0
            
            for i_row = 1:8
                for j_col = 1:12
                    
                    query_data = wildtype{i_row,j_col}{1};
                    log_query_data = log(query_data);
                    counts_list_query = histcounts(log_query_data,hist_bin_edge);
                    counts_total_query = size(query_data,1);
                    normalized_counts_query = counts_list_query ./ counts_total_query;
                    delta_freq = normalized_counts_query - normalized_counts_ref;
                    weight_on_peak = delta_freq;
                    weight_on_peak(weight_on_peak<0) = 0;
                    weight_on_peak(weight_on_peak>0 & [1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
                    
                    i_cutoff = find(weight_on_peak>0,1);
                    if i_cutoff
                        cutoff_wildtype(i_row,j_col) = hist_bin_edge(i_cutoff);
                    end
                    
                    ind_frac(i_row,j_col) = sum(weight_on_peak);
                    ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
                    basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
                    weight_off_peak = normalized_counts_query - weight_on_peak;
                    basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);
                    
                end
            end
            
            detect_limit = basal_frac(8,12);
            for i_row = 1:8
                for j_col = 1:12
                    if basal_frac(i_row,j_col) <= max(ind_frac(7,3), detect_limit)
                        mask_basal(i_row,j_col) = 0;
                        mask_induction(i_row,j_col) = 1;
                    elseif ind_frac(i_row,j_col) <= max(ind_frac(7,3), detect_limit)
                        mask_induction(i_row,j_col) = 0;
                        mask_basal(i_row,j_col) = 1;
                    end
                end
            end
            
    end         % end of switch
    
    % if isnan(basal or induced), use the other to fill up
    for i_row = 1:8
        for j_col = 1:12
            if isnan(basal_level(i_row,j_col))
                basal_level(i_row,j_col) = ind_level(i_row,j_col);
            end
            if isnan(ind_level(i_row,j_col))
                ind_level(i_row,j_col) = basal_level(i_row,j_col);
            end
        end
    end
    
    fitting_subset = {'all_data', '1r1c', '1r', '1c'};
    for i = 1:length(fitting_subset)
        fitting_type = fitting_subset{i};
        save_trait_mat( fitting_type, plate_name, basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc )
    end
    
end

%% sanity check on double gradient histogram
% dash lines show the cutoff

figure
set(gcf, 'position', [-1754 -196 1625 876]);
for i_row = 1:8
    for j_col = 1:12
        subplot(8,12,(i_row-1)*12+j_col)
        hist(log(mig1d{i_row,j_col}{1}),[0.05:0.1:9.95])
        hold on
        if ~isnan(cutoff_mig1d(i_row,j_col))
            addzeroline2('xpos',cutoff_mig1d(i_row,j_col));
        end
        xlim([0 10])
    end
end
export_fig(fullfile('../metaData/trait_extraction/', 'mig1d_double_gradient'))

figure
set(gcf, 'position', [-1754 -196 1625 876]);
for i_row = 1:8
    for j_col = 1:12
        subplot(8,12,(i_row-1)*12+j_col)
        hist(log(gal80d{i_row,j_col}{1}),[0.05:0.1:9.95])
        hold on
        if ~isnan(cutoff_gal80d(i_row,j_col))
            addzeroline2('xpos',cutoff_gal80d(i_row,j_col));
        end
        xlim([0 10])
    end
end
export_fig(fullfile('../metaData/trait_extraction/', 'gal80d_double_gradient'))

figure
set(gcf, 'position', [-1754 -196 1625 876]);
for i_row = 1:8
    for j_col = 1:12
        subplot(8,12,(i_row-1)*12+j_col)
        hist(log(wildtype{i_row,j_col}{1}),[0.05:0.1:9.95])
        hold on
        addzeroline2('xpos',cutoff_wildtype(i_row,j_col));
        xlim([0 10])
    end
end
export_fig(fullfile('../metaData/trait_extraction/', 'wildtype_double_gradient'))

