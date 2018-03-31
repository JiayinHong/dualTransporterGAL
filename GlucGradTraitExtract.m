%% load data

clear
load('../RawData/combined_filtered_data.mat')

%% analyze and extract trait data

% the level in the trait is log(yfp), use e as its base (natural logarithm)

gluc = [2 .*0.5 .^[0:9], 0, 0]';    %   use ' to convert the row array to column array
galc = [0, repmat(0.25,1,10), 2]';  %   repmat(A,M,N) is commonly used to produce an M-by-N
%   matrix filled with A's value and having A's CLASS.

strain_list = {'BC187', 'YJM978'};

for i_strain = 1:length(strain_list)
    strain_name = strain_list{i_strain};
    i_rows = find(strcmpi(strainData.name, strain_name));   %   strcmpi compare strings in a case insensitive way
    for i_rep = 1:length(i_rows)
        i_row = i_rows(i_rep);
        
        % prep trait
        basal_level = nan(12,1);
        ind_level = nan(12,1);
        modality = nan(12,1);
        basal_frac = nan(12,1);
        ind_frac = nan(12,1);
        mask_basal = nan(12,1);
        mask_induction = nan(12,1);
        
        switch strain_name
            case 'BC187'
                for i_col = 1:12
                    tmp = strainData{i_row, 'query'}{i_col};
                    tmp2 = log(tmp.yfp);
                    basal_level(i_col) = mean(tmp2);
                    ind_level(i_col) = basal_level(i_col);
                    modality(i_col) = 1;    %   BC187 is a unimodal strain
                    basal_frac(i_col) = 1;
                    ind_frac(i_col) = 0;
                    mask_basal(i_col) = 0.5;    %   the mask is hard-coded
                    mask_induction(i_col) = 0.5;
                end
                
            case 'YJM978'
                cutoff = [repmat(2.6,1,7), repmat(4,1,5)];  %   cutoff is hard-coded
                for i_col = 1:12
                    tmp = strainData{i_row, 'query'}{i_col};
                    tmp2 = log(tmp.yfp);
                    basal_level(i_col) = mean(tmp2(tmp2 < cutoff(i_col)));
                    ind_level(i_col) = mean(tmp2(tmp2 > cutoff(i_col)));
                    %   a little trick
                    basal_frac(i_col) = mean(tmp2 < cutoff(i_col));
                    ind_frac(i_col) = mean(tmp2 > cutoff(i_col));
                    %   tmp2 < cutoff(i_col) returns a logical array which
                    %   only contains number 0 or 1, so the mean of the
                    %   logical array is actually the fraction of the off peak
                end
                basal_level([11,12]) = ind_level([11,12]);
                ind_level([1:3]) = basal_level([1:3]);
                mask_basal = [...
                    0.5 0.5 0.5 1 ...
                    1   1   1   1 ...
                    1   0   0.5 0.5]';
                mask_induction = [...
                    0.5 0.5 0.5 0 ...
                    1   1   1   1 ...
                    1   1   0.5 0.5]';
                %   mask is hard-coded, (0.5,0.5) for a unimodal
                %   distribution; (1,0) or (0,1) for an uncertain
                %   distribution, in which only one peak is taken into
                %   account; and (1,1) for a bimodal distribution. Thus
                %   each spot we believe has a mask of 1
        end
        
        trait = table(basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc);
        
        save(...
            fullfile('../metaData/trait_extraction/', [strain_name, num2str(i_rep, '_rep_%02d')])...
            , 'trait');
        
        % plot yfp histogram
        [ yfp_hist_bin_edge, yfp_hist_bin_center ] = GetHistBin( 0, 10, 0.1 );
        figure
        for i_col = 1:12
            subplot(3,4,i_col)
            tmp = strainData{i_row, 'query'}{i_col};
            tmp2 = log(tmp.yfp);
            yfp_counts = histcounts(tmp2, yfp_hist_bin_edge);
            plot(yfp_hist_bin_center, yfp_counts, '-', 'LineWidth', 1.5)
            hold on
            addzeroline2('xpos', basal_level(i_col));
            addzeroline2('xpos', ind_level(i_col));
            xlim([0,10])
            xlabel('yfp')
            ylabel('counts')
        end
        export_fig(fullfile('../metaData/trait_extraction/', [strain_name, num2str(i_rep, '_rep_%02d')]))
    end
end