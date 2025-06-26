function [hmT] = getFirstHour(mT)
    % 22: active lever 22
    % 23: inactive lever 23
    % 95: head entries 95
    % 96: rewarded lever presses preceding head entries 96
    % 97: rewarded lever presses 97
    % 98: head entries following rewarded lever presses 98
    % 99: head entries (includes head entries following human-given reward + cue)
    % 17: infusion on 17
    % earned infusions count == rewarded leverpresses
    % total infusions count == infusion on
    
    copyVars = {'TagNumber', 'Session', 'sessionType', 'slideSession','Sex','TimeOfBehavior', 'LHbTarget', 'LHbAAV', 'Chamber', 'Acquire'};
    hmT = mT(:, copyVars);

    hourVars = {'ActiveLever', 'InactiveLever', 'HeadEntries', 'EarnedInfusions' ...
                'TotalInfusions', 'RewardedHeadEntries', 'RewardedLeverPresses'}; %SSnote: don't love the ambiguity of column name "RewardedLeverPresses"
    hourCodes = [22, 23, 95, 97, 17, 99, 97]; % codes for items in hourVars
    allEC = cell([height(mT), 1]);
    allET = cell([height(mT), 1]);
    for hv = 1:length(hourVars)
        dat = nan([height(mT), 1]);
        for fl = 1:height(mT)
            EC = mT.eventCode{fl};
            ET = mT.eventTime{fl};
            EC = EC(ET <= 3600);
            ET = ET(ET <= 3600); 

            if strcmp(hourVars{hv}, 'EarnedInfusions') | strcmp(hourVars{hv}, 'TotalInfusions')
                % calculate earned and total infusions 
                if mT.sessionType(fl) == categorical("Reinstatement")
                    dat(fl) = 0;
                elseif mT.TotalInfusions(fl) == 0
                    dat(fl) = 0;
                else
                    dat(fl) = length(find(EC==hourCodes(hv)));
                end
            else
                dat(fl) = length(find(EC==hourCodes(hv)));
            end
            
            % disp(length(find(EC==hourCodes(hv))))
            if strcmp(hourVars{hv}, 'InactiveLever') & (dat(fl) == 0)
                % SS hack for absent inactive lever event codes in some sessions...
                % HARDCODED TO READ RIGHT LEVER AS INACTIVE LEVER
                dat(fl) = length(find(EC==1));
            end
    
            if hv == 1
                allEC{fl} = EC;
                allET{fl} = ET;
            end
        end
        hmT = [hmT, table(dat)];
        hmT = renamevars(hmT, 'dat', hourVars{hv});
    end
    hmT = [hmT, table(allEC, allET)];
    hmT = renamevars(hmT, {'allEC', 'allET'}, {'eventCode', 'eventTime'});
    
    allLatency = cell([height(mT), 1]);
    Latency = nan([height(mT), 1]);
    doseHE = cell([height(mT), 1]);
    Intake = nan([height(mT), 1]);
    totalIntake = nan([height(mT), 1]);

    for fl = 1:height(mT)
        ET = hmT.eventTime{fl};
        EC = hmT.eventCode{fl};
        
        rewHE = ET(EC==99);
        doseHE{fl} = mT.doseHE{fl}(1:length(rewHE)); 

        if hmT.sessionType(fl) == 'Extinction'
            EC = hmT.eventCode{fl};
            ET = hmT.eventTime{fl};
            actLP = ET(EC==22);
            HE = ET(EC==95);
            seekHE = arrayfun(@(x) find(HE > x, 1, 'first'), actLP, 'UniformOutput', false);
            seekHE = HE(unique(cell2mat(seekHE(~cellfun(@isempty, seekHE)))));
            seekLP = arrayfun(@(x) find(actLP < x, 1, 'last'), seekHE, 'UniformOutput', false);
            seekLP = actLP(unique(cell2mat(seekLP(~cellfun(@isempty, seekLP)))));
            allLatency{fl} = seekHE-seekLP;
        else  
            allLatency{fl} = mT.allLatency{fl}(1:length(rewHE));
        end

        Latency(fl) = mean(allLatency{fl});
     
        if mT.TotalInfusions(fl) == 0
            totalIntake(fl) = 0;
            Intake(fl) = 0;
        else
            conc = mT.Concentration(fl);
            doseVol = mT.DoseVolume(fl);
            totalIntake(fl) = mT.TotalInfusions(fl) * doseVol * conc;
            Intake(fl) = mT.EarnedInfusions(fl) * doseVol * conc;
        end
    end
    hmT = [hmT, table(allLatency, Latency, doseHE, Intake, totalIntake)];
end