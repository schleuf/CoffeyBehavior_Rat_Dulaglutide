function [PCA] = PCA_analysis(ivZT, pcaGroups, sub_dir, saveTabs, tabs_savepath, z_suff)
    PCA = struct; 
    for pg = 1:length(pcaGroups)
        use_inds = ones([height(ivZT), 1]);
        if strcmp(pcaGroups{pg}{1}{1}, 'all')
            suff_str = [z_suff, '_all']; 
        else
            suff_str = z_suff;
            for cat = 1:length(pcaGroups{pg})
                suff_str = [suff_str, '_', pcaGroups{pg}{cat}{2}];
                use_inds = use_inds & (ivZT.(pcaGroups{pg}{cat}{1})==(pcaGroups{pg}{cat}{2}));
            end
        end
        if length(find(use_inds)) > 3
            this_field = suff_str(2:end);
            PCA.(this_field) = struct;
            PCA.(this_field).ivZT_inds = find(use_inds);   
            prednames = ivZT.Properties.VariableNames;
            prednames = prednames(~ismember(prednames, {'ID', 'Strain', 'Sex', 'Acquire', 'Severity', 'Class'}));
            [coeff,score,latent] = pca(ivZT{PCA.(this_field).ivZT_inds, prednames});
            PCA.(this_field).coeff = coeff;
            PCA.(this_field).score = score;
            PCA.(this_field).latent = latent;
            PCA.(this_field).prednames = prednames;
        else
            disp(['not enough data available for PCA: ', suff_str(2:end)])
        end
    end
    if saveTabs
        save([sub_dir, tabs_savepath, 'IS_metric_PCA', char(z_suff), '.mat'], 'PCA');
    end
end