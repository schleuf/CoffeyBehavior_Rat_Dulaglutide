function [varTable] = rawVariableExtractor(varTable, eventCode, eventTime)
% rawVariableExtractor takes the raw event times for mouse oral SA and
% calculates variables of interest and adds them to the varTable
%
% INPUTS: varTable, eventCode, eventTime - Direct outputs of "importMouseOralSA"
%
% OUTPUTS: varTable (Updated with variables calculated or cleaned below)
%
% Variables added to or modified from varTable:
%    HeadEntries: # head entries (filtered to remove headentries <2s apart)
%    RewardedHeadEntries: # head entries occurring after a lever press 
%    RewardedLeverPresses: # active lever presses with a head entry before the next activelever press
%    DoseHE: # of infusions per rewarded head entry
%    allLatency: latencies from rewarded lever presses to first head entry (unless no head entry occurs before next rewarded press)
%    Latency: Avg. time from each rewarded lever press to first head entry
%    eventCode: original eventCode appended with new event codes (see below)
%    eventTime: original eventTime appended with new event times (see below)
% 
% events added to eventCode and eventTime:
%    95 time-filtered head entries
%    96 rewarded lever presses preceding head entries 
%    97 rewarded lever presses (includes lever presses not followed by a head entry before the next reward)
%    98 rewarded head entries following active lever presses
%    99 rewarded head entries (includes head entries following human-given reward + cue)
    
    showWarnings = false;

    % Change Subject to TagNumber to Match the Key
    varTable.Properties.VariableNames{'Subject'} = 'TagNumber';
    
    % KRCnote: UPDATE currently doesn't work with PR data -- array length of eventCode >
    % eventTime, so I just made a small work around so I can see ~how the
    % animals are doing during PR
    ldiff=length(eventCode)-length(eventTime);
    if ldiff > 0
        eventCode = eventCode(1:end-ldiff);
    end

    HE_timeFilt = 2;
    % filter head entries, remove entries that happen within  HE_timeFilt(s) of each other
    time_prefiltHE=eventTime(eventCode==6);
    remove_inds = logical([0; diff(time_prefiltHE) < HE_timeFilt]);
    time_HE = time_prefiltHE;
    time_HE(remove_inds) = [];
    varTable.HeadEntries = length(time_HE); 
    
    % filter for rewarded head entries (the first occurring after a lever press)
    time_cue = eventTime(eventCode==13); 
    nextHE = arrayfun(@(x) time_HE(find(time_HE > x, 1, 'first')), time_cue, 'UniformOutput', false);
    nextHE = cell2mat(nextHE(~cellfun(@isempty, nextHE))); % Remove empty cells and convert to numeric
    time_rewHE = unique(nextHE, 'stable'); % Ensure uniqueness
    varTable.rewardedHeadEntries = length(time_rewHE);

    time_rewLP = eventTime(eventCode == 3 | eventCode == 4);
    % filter out rewarded lever presses not followed by head entries and
    % rewarded head entries not preceded by active lever presses (freebies)
    time_rewLP_preceding_HE = arrayfun(@(x) time_rewLP(find(time_rewLP < x, 1, 'last')), time_rewHE, 'UniformOutput', false);
    time_rewLP_preceding_HE = cell2mat(time_rewLP_preceding_HE(~cellfun(@isempty, time_rewLP_preceding_HE))); 
    time_rewLP_preceding_HE = unique(time_rewLP_preceding_HE, 'stable');
    
    time_HE_following_rewLP = arrayfun(@(x) time_rewHE(find(time_rewHE > x, 1, 'first')), time_rewLP, 'UniformOutput', false);
    time_HE_following_rewLP = cell2mat(time_HE_following_rewLP(~cellfun(@isempty, time_HE_following_rewLP)));
    time_HE_following_rewLP = unique(time_HE_following_rewLP, 'stable');
    
    % filter for the lever presses that precede the rewarded head entries and get dose per head entry
    doseHE = arrayfun(@(x) sum(time_rewLP < x), time_HE_following_rewLP); % Compute the number of rewLP timestamps before each rewHE timestamp
    if ~isempty(doseHE)
        doseHE = [doseHE(1); diff(doseHE)]; % Compute the difference between successive elements to get the count per interval
    end
    varTable.doseHE = {doseHE};

    % calculate earned and total infusions 
    if varTable.Session == 26 % SSnote: maybe this get's moved to createNewMasterTable too
        varTable.EarnedInfusions = 0;
    elseif varTable.TotalInfusions == 0
        varTable.EarnedInfusions = 0;
    else
        varTable.EarnedInfusions = length(time_rewLP);
    end

    % calculate all and mean latency
    varTable.allLatency = {time_HE_following_rewLP - time_rewLP_preceding_HE};
    trimLatency=varTable.allLatency{1};
    varTable.Latency = mean(trimLatency(trimLatency<360)); % Setting Max Latency to 360s 
    varTable.Latency(isempty(varTable.allLatency{1})) = NaN;

    % if (varTable.Session > 15) && (varTable.Session <=26)
    %     disp('bleh');
    % end
    
    % get number of active and inactive lever presses that occur during ITI
    % SSnote: these are encoded as 20 and 21, don't need to calc them 
    varTable.itiActiveLever= varTable.ActiveLever - varTable.EarnedInfusions;
    inLP=eventTime(eventCode==23);% 23  = Inactive Lev press
    % % SS hack for absent inactive lever event codes in some sessions...
    % if isempty(inLP)
    %     inLP = eventTime(eventCode==1); % HARDCODED TO READ RIGHT LEVER AS INACTIVE LEVER
    %     if showWarnings
    %         if ~isempty(inLP)
    %             disp(['logged right lever presses as inactive lever presses for ID:', char(varTable.TagNumber), ' on ', char(varTable.Date)])
    %         else
    %             disp(['no inactive lever presses found for ID:', char(varTable.TagNumber), ' on ', char(varTable.Date)])
    %         end
    %     end
    % end
    inITI=[]; % Initialize in case there are no inLPs during ITI   
    for j = 1:height(time_rewLP)
        inITI(j,1)=sum(inLP>time_rewLP(j) & inLP<(time_rewLP(j)+10)); % Which inactive LPs occur in the 10s following a reward (the ITI)
    end
    varTable.itiInactiveLever=sum(inITI);
    
    % append eventCode and eventTime for new variables 
    eventCode=[eventCode; ...
               repmat(95, [height(time_HE), 1]); ...
               repmat(96,[height(time_rewLP),1]); ... 
               repmat(97,[height(time_rewLP_preceding_HE),1]); ...
               repmat(98,[height(time_rewHE),1]); ...
               repmat(99,[height(time_HE_following_rewLP),1])];
    
    eventTime = [eventTime; ...
                 time_HE; ...
                 time_rewLP; ... % lever presses preceding earned infusion
                 time_rewLP_preceding_HE; ...
                 time_rewHE; ... % head entries following cue light
                 time_HE_following_rewLP];

    varTable.eventCode={eventCode};
    varTable.eventTime={eventTime};
end