function [ivT] = IS_processes(mT, dex, runType, corrGroups, violSubsets, violGroups, violLabels, pcaGroups, sub_dir, saveTabs, tabs_savepath, groupOralFentOutput_figs, groupOralFentOutput_savepath, figsave_type)
    
    % get IS metrics
    [ivT] = GetMetrics(mT(dex.all, :));
    if saveTabs
        save([sub_dir, tabs_savepath, 'IS_Metrics', '.mat'], 'ivT');
    end

    % run Z scores separately for animals that lack ER sessions & associated metrics
    if ~any(ismember(runType, 'ER'))
        zGroups = {ivT.Acquire=='Acquire'}; %{ones([height(ivT), 1])};
        includeER = false;
        z_suff = {'_noER'};
    elseif (length(runType) > 1) & any(ismember(runType, 'ER'))
        ER_IDs =  unique(mT.TagNumber(dex.ER));
        ivT_ER_ind = ismember(unique(mT.TagNumber(dex.all)), ER_IDs);
        zGroups = {ivT.Acquire=='Acquire', ivT_ER_ind & ivT.Acquire=='Acquire'}; %{ones([height(ivT), 1]), ivT_ER_ind};
        includeER = [false, true];
        z_suff = {'_noER', '_withER'};
    else % only 'ER' in runType
        zGroups = {ivT.Acquire=='Acquire'}; %{ones([height(ivT), 1])};
        includeER = true;
        z_suff = {'_withER'};
    end

    for zg = 1:length(zGroups)
        [ivZT, removed_ind] = SeverityScore(ivT(logical(zGroups{zg}),:), includeER(zg));
        if saveTabs
            save([sub_dir, tabs_savepath, 'IS_Zscores', char(z_suff{zg}), '.mat'], 'ivZT');
        end
        [correlations] = GetCorr(ivZT, z_suff{zg}, corrGroups, sub_dir, groupOralFentOutput_figs, groupOralFentOutput_savepath, figsave_type);
        if groupOralFentOutput_figs     
            for vg = 1:length(violSubsets)
                thistab = ivT(logical(zGroups{zg}), :);
                thistab = thistab(~removed_ind, :);
                thistab.Severity = ivZT.Severity;
                subset = ones([height(thistab), 1]);
                if ~strcmp(violSubsets{vg}{1}, 'all')
                    subset = subset & (thistab.(violSubsets{vg}{1}) == violSubsets{vg}{2});
                end
                % if ~isempty(find(subset, 1)) && length(unique(thistab(subset,:).(violGroups{vg}))) > 1
                %     ViolinFig(thistab(logical(subset), :), violGroups{vg}, [z_suff{zg}, '_', violLabels{vg}], includeER(zg), sub_dir, groupOralFentOutput_savepath, figsave_type)
                % else
                try
                    ViolinFig(thistab(logical(subset), :), violGroups{vg}, [z_suff{zg}, '_', violLabels{vg}], includeER(zg), sub_dir, groupOralFentOutput_savepath, figsave_type, violLabels{vg})
                catch
                    disp(['no data available for violin plot: ', z_suff{zg}, '_', violLabels{vg}])
                end
            end
            PCA = PCA_analysis(ivZT, pcaGroups, sub_dir, saveTabs, tabs_savepath, z_suff{zg});
            if groupOralFentOutput_figs 
                PCA_fig(ivZT, PCA, sub_dir, groupOralFentOutput_savepath, z_suff{zg}, figsave_type);
            end
        end
    end
end