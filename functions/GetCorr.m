function [correlations] = GetCorr(ivZT, z_suff, corrGroups, sub_dir, groupOralFentOutput_figs, groupOralFentOutput_savepath, figsave_type)
    correlations = struct; 
    for cg = 1:length(corrGroups)
        use_inds = ones([height(ivZT), 1]);
        if strcmp(corrGroups{cg}{1}{1}, 'all')
            suff_str = [z_suff, '_all']; 
        else
            suff_str = z_suff;
            for cat = 1:length(corrGroups{cg})
                suff_str = [suff_str, '_', corrGroups{cg}{cat}{2}];
                use_inds = use_inds & (ivZT.(corrGroups{cg}{cat}{1})==(corrGroups{cg}{cat}{2}));
            end
        end
        if ~isempty(find(use_inds,1))
            correlations.(suff_str(2:end)) = struct;
            correlations.(suff_str(2:end)).ivZT_inds = find(use_inds);   
            prednames = ivZT.Properties.VariableNames;
            prednames = prednames(~ismember(prednames, {'ID', 'Strain', 'Sex', 'Acquire', 'Severity', 'Class'}));
            correlations.(suff_str(2:end)).ct = corr(ivZT{find(use_inds),prednames},Type='Pearson');
            correlations.(suff_str(2:end)).prednames = prednames;
            if groupOralFentOutput_figs
                CorrFig(correlations.(suff_str(2:end)).ct, prednames, sub_dir, groupOralFentOutput_savepath, suff_str, figsave_type)
            end
        else

            disp(['no data available for correlations in this group: ' suff_str(2:end)])
    end
end