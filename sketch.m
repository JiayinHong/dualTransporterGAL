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