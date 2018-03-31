%% load data

clear
load('../RawData/segmented-GAL234pr-YFP-data.mat')

% the dimensions of 'strainData': first dimension - strain type,
% 1 for GAL2pr, 2 for GAL3pr, 3 for GAL4pr; second dimension - glucose gradient;
% third dimension - galactose gradient;

% GAL2promoter fused with YFP
gal2pr = cell(8,12);
gal2max = nan(8,12);
gal2min = nan(8,12);
for i_row = 1:8
	for j_col = 1:12
		tmp = strainData(1,i_row,j_col).yfp;
		gal2pr{i_row,j_col} = tmp;
        gal2max(i_row,j_col) = max(tmp);
        gal2min(i_row,j_col) = min(tmp);
	end
end

% GAL3promoter fused with YFP
gal3pr = cell(8,12);
gal3max = nan(8,12);
gal3min = nan(8,12);
for i_row = 1:8
	for j_col = 1:12
		tmp = strainData(2,i_row,j_col).yfp;
		gal3pr{i_row,j_col} = tmp;
        gal3max(i_row,j_col) = max(tmp);
        gal3min(i_row,j_col) = min(tmp);
	end
end

% GAL4promoter fused with YFP
gal4pr = cell(8,12);
gal4max = nan(8,12);
gal4min = nan(8,12);
for i_row = 1:8
	for j_col = 1:12
		tmp = strainData(3,i_row,j_col).yfp;
		gal4pr{i_row,j_col} = tmp;
        gal4max(i_row,j_col) = max(tmp);
        gal4min(i_row,j_col) = min(tmp);
	end
end

% prepare to plot histogram
galLabel = {'None','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2'};
gluLabel = {'0','-1','-2','-3','-4','-5','-6','None'};
% colormap = parula;
colormap = jet;

% mean metric
m0 = 1.3;
mmax = 4;

%% make the trait table for double gradient condition, GAL2pro, GAL3pro, and GAL4pro

for plate_type = {'GAL2pr', 'GAL3pr', 'GAL4pr'}
	plate_type = plate_type{1}; 	% convert cell to string

	% prep trait
	basal_level = nan(8,12);
	ind_level = nan(8,12);
	modality = nan(8,12);
	basal_frac = nan(8,12);
	ind_frac = nan(8,12);
	mask_basal = 0.5 .* ones(8,12);
	mask_induction = 0.5 .* ones(8,12);

	gluc_gradient = [2 .^[0:-1:-6], 0]';
	galc_gradient = [0, 2 .^[-8:1:2]];
	gluc = gluc_gradient * ones(1,12);
	galc = ones(8,1) * galc_gradient;

% 	[hist_bin_edge, hist_bin_center] = GetHistBin(0.05, 9.95, 0.1);
    [hist_bin_edge, hist_bin_center] = GetHistBin(0.05, 8.05, 0.1);


	switch plate_type

	case 'GAL2pr'

		cutoff_gal2pr = nan(8,12);
		ref_data = gal2pr{1,1};	% fully OFF
		log_ref_data = log(ref_data);	% natural logarithm
		counts_total_ref = numel(ref_data);	% the total counts of ref data
		% get a list of counts in each bin
		counts_list_ref = histcounts(log_ref_data, hist_bin_edge);
		normalized_counts_ref = counts_list_ref ./ counts_total_ref;
		threshold = mean(log_ref_data);
		% ignore if there are more normalized counts of query than ref data
		% beneath the threshold

		[hf ha] = gridplot(8,12,70,70,'gap',2);

		for i_row = 1:8
			for j_col = 1:12

				query_data = gal2pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				% calculate the difference between query and fully OFF, 
				% regarded as the ON fraction of query data
				delta_freq = normalized_counts_query - normalized_counts_ref;
				% get the weight of ON peak in each bin
				weight_on_peak = delta_freq;
				weight_on_peak(weight_on_peak<0) = 0;
				weight_on_peak(weight_on_peak>0 & ...
					[1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
				% 0.1 is the bin size
				i_cutoff = find(weight_on_peak>0, 1);
				% find in which bin that the first ON fraction lies in
				if i_cutoff
					cutoff_gal2pr(i_row,j_col) = hist_bin_edge(i_cutoff);
				end
				% get the cutoff line position

				ind_frac(i_row,j_col) = sum(weight_on_peak);
				ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
				basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
				weight_off_peak = normalized_counts_query - weight_on_peak;
				basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);

				% plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
% 				ax.Color = colormap(min(length(colormap), floor((ind_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
% 				ax.XTickLabel = [];
				ax.YTickLabel = [];

				if i_row==8		% the last row need label
					xlabel(galLabel{j_col}, 'FontSize', 20);
				end
				if j_col==1		% the first col need label
					h = ylabel(gluLabel{i_row}, 'FontSize', 20);
					h.Rotation = 0;
					h.HorizontalAlignment = 'right';
					h.VerticalAlignment = 'middle';
					h.Margin = 10;
                end
                

			end
        end
        
        [ax,h] = suplabel('GAL2pro-YFP', 't');
        h.FontSize = 20;

		detect_limit = basal_frac(8,12);
		% the lower right corner should be fully ON, where basal_frac should be zero

		for i_row = 1:8
			for j_col = 1:12
				if basal_frac(i_row,j_col) <= detect_limit
					mask_basal(i_row,j_col) = 0;
					% mask_induction(i_row,j_col) = 1;
				elseif ind_frac(i_row,j_col) <= detect_limit
					mask_induction(i_row,j_col) = 0;
					% mask_basal(i_row,j_col) = 1;
				end
			end
		end


	case 'GAL3pr'

		cutoff_gal3pr = nan(8,12);
		ref_data = gal3pr{1,1};	% fully OFF
		log_ref_data = log(ref_data);	% natural logarithm
		counts_total_ref = numel(ref_data);	% the total counts of ref data
		% get a list of counts in each bin
		counts_list_ref = histcounts(log_ref_data, hist_bin_edge);
		normalized_counts_ref = counts_list_ref ./ counts_total_ref;
		threshold = mean(log_ref_data);
		% ignore if there are more normalized counts of query than ref data
		% beneath the threshold

        [hf ha] = gridplot(8,12,70,70,'gap',2);
        
		for i_row = 1:8
			for j_col = 1:12

				query_data = gal3pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				% calculate the difference between query and fully OFF, 
				% regarded as the ON fraction of query data
				delta_freq = normalized_counts_query - normalized_counts_ref;
				% get the weight of ON peak in each bin
				weight_on_peak = delta_freq;
				weight_on_peak(weight_on_peak<0) = 0;
				weight_on_peak(weight_on_peak>0 & ...
					[1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
				% 0.1 is the bin size
				i_cutoff = find(weight_on_peak>0, 1);
				% find in which bin that the first ON fraction lies in
				if i_cutoff
					cutoff_gal3pr(i_row,j_col) = hist_bin_edge(i_cutoff);
				end
				% get the cutoff line position

				ind_frac(i_row,j_col) = sum(weight_on_peak);
				ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
				basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
				weight_off_peak = normalized_counts_query - weight_on_peak;
				basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);
                
                % plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
% 				ax.Color = colormap(min(length(colormap), floor((ind_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
% 				ax.XTickLabel = [];
				ax.YTickLabel = [];

				if i_row==8		% the last row need label
					xlabel(galLabel{j_col}, 'FontSize', 20);
				end
				if j_col==1		% the first col need label
					h = ylabel(gluLabel{i_row}, 'FontSize', 20);
					h.Rotation = 0;
					h.HorizontalAlignment = 'right';
					h.VerticalAlignment = 'middle';
					h.Margin = 10;
				end


			end
        end

        [ax,h] = suplabel('GAL3pro-YFP', 't');
        h.FontSize = 20;
        
		detect_limit = basal_frac(8,12);
		% the lower right corner should be fully ON, where basal_frac should be zero

		for i_row = 1:8
			for j_col = 1:12
				if basal_frac(i_row,j_col) <= detect_limit
					mask_basal(i_row,j_col) = 0;
					% mask_induction(i_row,j_col) = 1;
				elseif ind_frac(i_row,j_col) <= detect_limit
					mask_induction(i_row,j_col) = 0;
					% mask_basal(i_row,j_col) = 1;
				end
			end
		end

	case 'GAL4pr'

		cutoff_gal4pr = nan(8,12);
		ref_data = gal4pr{1,1};	% fully OFF
		log_ref_data = log(ref_data);	% natural logarithm
		counts_total_ref = numel(ref_data);	% the total counts of ref data
		% get a list of counts in each bin
		counts_list_ref = histcounts(log_ref_data, hist_bin_edge);
		normalized_counts_ref = counts_list_ref ./ counts_total_ref;
		threshold = mean(log_ref_data);
		% ignore if there are more normalized counts of query than ref data
		% beneath the threshold
        
        [hf ha] = gridplot(8,12,70,70,'gap',2);

		for i_row = 1:8
			for j_col = 1:12

				query_data = gal4pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				% calculate the difference between query and fully OFF, 
				% regarded as the ON fraction of query data
				delta_freq = normalized_counts_query - normalized_counts_ref;
				% get the weight of ON peak in each bin
				weight_on_peak = delta_freq;
				weight_on_peak(weight_on_peak<0) = 0;
				weight_on_peak(weight_on_peak>0 & ...
					[1:length(weight_on_peak)] <= ceil(threshold/0.1)) = 0;
				% 0.1 is the bin size
				i_cutoff = find(weight_on_peak>0, 1);
				% find in which bin that the first ON fraction lies in
				if i_cutoff
					cutoff_gal4pr(i_row,j_col) = hist_bin_edge(i_cutoff);
				end
				% get the cutoff line position

				ind_frac(i_row,j_col) = sum(weight_on_peak);
				ind_level(i_row,j_col) = sum(hist_bin_center .* weight_on_peak) / sum(weight_on_peak);
				basal_frac(i_row,j_col) = 1-ind_frac(i_row,j_col);
				weight_off_peak = normalized_counts_query - weight_on_peak;
				basal_level(i_row,j_col) = sum(hist_bin_center .* weight_off_peak) / sum(weight_off_peak);
                
                % plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
% 				ax.Color = colormap(min(length(colormap), floor((ind_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
% 				ax.XTickLabel = [];
				ax.YTickLabel = [];

				if i_row==8		% the last row need label
					xlabel(galLabel{j_col}, 'FontSize', 20);
				end
				if j_col==1		% the first col need label
					h = ylabel(gluLabel{i_row}, 'FontSize', 20);
					h.Rotation = 0;
					h.HorizontalAlignment = 'right';
					h.VerticalAlignment = 'middle';
					h.Margin = 10;
				end


			end
        end

        [ax,h] = suplabel('GAL4pro-YFP', 't');
        h.FontSize = 20;
        
		detect_limit = basal_frac(8,12);
		% the lower right corner should be fully ON, where basal_frac should be zero

		for i_row = 1:8
			for j_col = 1:12
				if basal_frac(i_row,j_col) <= detect_limit
					mask_basal(i_row,j_col) = 0;
					% mask_induction(i_row,j_col) = 1;
				elseif ind_frac(i_row,j_col) <= detect_limit
					mask_induction(i_row,j_col) = 0;
					% mask_basal(i_row,j_col) = 1;
				end
			end
		end

	end		% end of switch

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

% 	fitting_subset = {'all_data', '1r1c', '1r', '1c'};
    fitting_subset = {'all_data'};
    for i = 1:length(fitting_subset)
        fitting_type = fitting_subset{i};
        save_trait_mat( fitting_type, plate_type, basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc )
    end
    
end
						







