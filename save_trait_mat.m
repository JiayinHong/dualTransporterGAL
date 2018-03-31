function save_trait_mat( fitting_type, plate_name, basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc )
%   This function is called by 'DoubleGradTraitExtract', given a
%   fitting_type, generate corresponding mask values, and save into a
%   trait.mat file
switch fitting_type
    
    case 'all_data'
        mask_basal(8,1:4) = 0;
        mask_induction(8,1:4) = 0;
        
    case '1r1c'
        mask_basal([1:3,5:8],[1:8,10:12]) = 0;
        mask_induction([1:3,5:8],[1:8,10:12]) = 0;
        
    case '1r'
        mask_basal([1:3,5:8],:) = 0;
        mask_induction([1:3,5:8],:) = 0;
        
    case '1c'
        mask_basal(:,[1:8,10:12]) = 0;
        mask_induction(:,[1:8,10:12]) = 0;
        
end

trans_table_trait_extract;  % convert all the matrix to one column

trait = table(basal_level, ind_level, modality, basal_frac, ind_frac, mask_basal, mask_induction, gluc, galc);
save(fullfile('../metaData/trait_extraction', [plate_name,'_',fitting_type]), 'trait');

end