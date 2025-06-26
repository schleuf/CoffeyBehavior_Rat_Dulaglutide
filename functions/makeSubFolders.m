function [new_dirs] = makeSubFolders(allfig_savefolder, runNum, runType, toMake, excludeData, firstHour)
    
    if ~exist(allfig_savefolder, 'dir')
        mkdir(allfig_savefolder)
        disp(['Folder created: ', allfig_savefolder])
    end

    runTypeStr = 'Exp_';
    if length(runType) > 1
        runTypeStr = [runTypeStr, 'all'];
    else
        runTypeStr = [runTypeStr, char(string(runType))];
    end
    
    runNumStr = 'Run_';
    if length(runNum) > 1
        runNumStr = [runNumStr, 'all'];
    else
        runNumStr = [runNumStr, char(string(runNum))];
    end
    
    sub_dir = [runTypeStr, '_', runNumStr];
    if excludeData
        sub_dir = [sub_dir, '_exclusions'];
    end
    sub_dir = [allfig_savefolder, sub_dir, '\'];
    new_dirs = {sub_dir};

    if firstHour
        fH_sub_dir = [sub_dir(1:length(sub_dir)-1), '_firstHour', '\'];
        new_dirs = {new_dirs{1}, fH_sub_dir};
    end
    
    for tm = 1:length(toMake)
        if ~exist([sub_dir, toMake{tm}], 'dir')
            mkdir([sub_dir, toMake{tm}])
            disp(['Folder created: ', sub_dir, toMake{tm}])
        end
        if firstHour && ~exist([fH_sub_dir, toMake{tm}], 'dir')
            mkdir([fH_sub_dir, toMake{tm}])
            disp(['Folder created: ', fH_sub_dir, toMake{tm}])
        end
    end
end
