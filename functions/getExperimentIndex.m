function [dex] = getExperimentIndex(mT, runNum, runType)
    if runNum == 'all'
        runNum_inds = ones([height(mT), 1]);
    else
        nums = strsplit(char(runNum),'_');
        runNum_inds = zeros([height(mT),1]);
        for n = 1:length(nums)
            runNum_inds = runNum_inds | (mT.Run == str2double(nums{n}));
        end
    end
    
    % Create indexing structure for 'ER', 'BE', and 'SA' data of the desired run #(s)
    dex = struct;
    dex.all = [];
    for e = 1:length(runType)
        if runType(e) == 'SA'
            dex.SA = find(((mT.sessionType == 'Training') | (mT.sessionType == 'PreTraining')) & runNum_inds);
        else
            dex.(string(runType(e))) = find((mT.Experiment == runType(e)) & runNum_inds);
        end
        dex.all = union(dex.all, dex.(string(runType(e))));
    end
end
