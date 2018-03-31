%% export GAL1 fcs data
load('../RawData/mig1d_gal80d_wt_data.mat')
wildtype = cell(8,12);
for i_row = 8:-1:1
    for j_col = 1:12
        tmp = {alldata(5,i_row,j_col,2).yfp};
        wildtype{9-i_row, j_col} = tmp;
        filename = sprintf('well %s-%s.csv', num2str(9-i_row),num2str(j_col));
        csvwrite(fullfile('../GAL1csv-files/', filename),tmp)
    end
end

%% export GAL3 fcs data
load('../RawData/segmented-GAL234pr-YFP-data.mat')
gal3pr = cell(8,12);
for i_row = 1:8
	for j_col = 1:12
		tmp = strainData(2,i_row,j_col).yfp;
		gal3pr{i_row,j_col} = tmp;
        filename = sprintf('well %s-%s.csv', num2str(i_row),num2str(j_col));
        csvwrite(fullfile('../GAL3csv-files/', filename),tmp)
	end
end