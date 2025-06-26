function [Acquire] = getAcquire(mT, acquisition_thresh, acquisition_testPeriod, pAcq)
    IDs=unique(mT.TagNumber);
    mT = sortrows(mT,'TagNumber','ascend');
    Acq = nan(size(IDs));
    Acquire = categorical(nan([height(mT), 1]));
    
    sessionType = acquisition_testPeriod{1};
    sessionRange = acquisition_testPeriod{2};

    testSessions = unique(mT.Session(mT.sessionType == sessionType));
    if ~strcmp(sessionRange, 'all')
        if length(acquisition_testPeriod) > 2
            numGrab = acquisition_testPeriod{3};
        else
            numGrab = 1;
        end
        if strcmp(sessionRange, 'first')
            testSessions = testSessions(1:numGrab);
        elseif strcmp(sessionRange, 'last')
            testSessions = testSessions(end-numGrab+1:end);
        else
            disp('you used an input for acquisition test_period other than "all", "first", or "last"...defaulting to "all"')
        end
    end


    for id=1:height(IDs)
        % SSnote: this should be using sessionType instead of the Session#
        % idx = mT.TagNumber==IDs(id) & mT.Session>9 & mT.Session <16;
        idx = mT.TagNumber==IDs(id) & ismember(mT.Session, testSessions);
        Acq(id,1) = mean(mT.EarnedInfusions(idx)) > acquisition_thresh;
        acqDist(id,1) = mean(mT.EarnedInfusions(idx));
        
        disp(IDs(id))
        if ~Acq(id,1)
            
            disp(acqDist(id,1))
            if isnan(acqDist(id,1))
                disp(mT.EarnedInfusions(idx))
            end
            
        end
        disp(' ')

        if Acq(id) == 0 % && sum(idx) ~= 0
            status = 'NonAcquire';
        else
            status = 'Acquire';
        end
        tmp=repmat(categorical({status}), sum(mT.TagNumber == IDs(id)), 1);
        Acquire(mT.TagNumber == IDs(id)) = tmp;
    end
    if pAcq
       figure('Position',[1 300 275 250],'color','w'); histogram(acqDist,50);
       box off; hold on;
       plot([acquisition_thresh acquisition_thresh],[0 15],'--r')
       ylabel('Number of Animals');
       xlabel('Rewards Earned (Training Avg)');
       set(gca,'LineWidth',1.5,'TickDir','out','FontSize',12);
    end
end