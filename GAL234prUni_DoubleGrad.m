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
m0 = 1.28;
mmax = 3.8;

%% make the trait table for double gradient condition, GAL2pro, GAL3pro, and GAL4pro

for plate_type = {'GAL2pr', 'GAL3pr', 'GAL4pr'}
	plate_type = plate_type{1}; 	% convert cell to string

	% prep trait
	basal_level = nan(8,12);
	ind_level = nan(8,12);
	modality = nan(8,12);
	basal_frac = nan(8,12);
	ind_frac = nan(8,12);
	mask_basal = ones(8,12);
	mask_induction = ones(8,12);

	gluc_gradient = [2 .^[0:-1:-6], 0]';
	galc_gradient = [0, 2 .^[-8:1:2]];
	gluc = gluc_gradient * ones(1,12);
	galc = ones(8,1) * galc_gradient;

% 	[hist_bin_edge, hist_bin_center] = GetHistBin(0.05, 9.95, 0.1);
    [hist_bin_edge, hist_bin_center] = GetHistBin(0.05, 8.05, 0.1);


	switch plate_type

	case 'GAL2pr'

		[hf ha] = gridplot(8,12,70,70,'gap',2);

		for i_row = 1:8
			for j_col = 1:12

				query_data = gal2pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				basal_level(i_row,j_col) = mean(log_query_data);

				% plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
				ax.Color = colormap(min(length(colormap), floor((basal_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
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


	case 'GAL3pr'

        [hf ha] = gridplot(8,12,70,70,'gap',2);
        
		for i_row = 1:8
			for j_col = 1:12

				query_data = gal3pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				basal_level(i_row,j_col) = mean(log_query_data);
                
                % plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
				ax.Color = colormap(min(length(colormap), floor((basal_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
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
        

	case 'GAL4pr'
        
        [hf ha] = gridplot(8,12,70,70,'gap',2);

		for i_row = 1:8
			for j_col = 1:12

				query_data = gal4pr{i_row, j_col};
				log_query_data = log(query_data); 	% natural logarithm
				counts_total_query = numel(query_data);
				counts_list_query = histcounts(log_query_data, hist_bin_edge);
				normalized_counts_query = counts_list_query ./ counts_total_query;

				basal_level(i_row,j_col) = mean(log_query_data);
                
                % plot
				axes(ha(sub2ind([12 8],j_col,i_row)));
				plot(hist_bin_center, normalized_counts_query, 'k-', 'LineWidth', 2);
                xlim([0 8.1])
				hold all
				
				drawnow
				ax = gca;
				ax.Color = colormap(min(length(colormap), floor((basal_level(i_row,j_col)-m0)/mmax*length(colormap))+1),:);
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

	end		% end of switch


% 	fitting_subset = {'all_data', '1r1c', '1r', '1c'};
    fitting_subset = {'all_data'};
    for i = 1:length(fitting_subset)
        fitting_type = fitting_subset{i};
        save_trait_mat( fitting_type, plate_type, basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc )
    end
    
end
						







