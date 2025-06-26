function [ivZT, remove_rows] = SeverityScore(ivT, includeER)
    
    ivZT = ivT(:, {'ID', 'Sex', 'Strain', 'Acquire'});

    % all metrics are measured from training epoch. remove animals that
    % haven't done the training epoch yet. 
    remove_rows = isnan(ivT.Intake);
    if ~isempty(find(remove_rows, 1))
        disp('removed animals from non-ER PCA: ')
        disp(ivT.ID(remove_rows)');
    end
        
    % get Z scores for non-ER IS metrics
    [int, seek, asso, esc] = deal(nan([height(ivT), 1]));
    int(~remove_rows) = zscore(ivT.Intake(~remove_rows));
    seek(~remove_rows) = zscore(ivT.Seeking(~remove_rows));
    asso(~remove_rows) = zscore(ivT.Association(~remove_rows));
    esc(~remove_rows) = zscore(ivT.Escalation(~remove_rows));
    ivZT.Intake = int;
    ivZT.Seeking = seek;
    ivZT.Association = asso;
    ivZT.Escalation = esc;

    if includeER
        % remove animals that haven't done the escalation or relapse epochs yet. 
        remove_rows = remove_rows | isnan(ivT.Escalation) | isnan(ivT.Relapse);
        if ~isempty(find(remove_rows, 1))
            disp('removed animals from ER PCA: ')
            disp(ivT.ID(find(isnan(ivT.Escalation) | isnan(ivT.Relapse)))');
        end
        [esc, rel, rec] = deal(nan([height(ivT), 1]));
        % get Z scores for ER IS metrics
        esc(~remove_rows) = zscore(ivT.Extinction(~remove_rows));
        rel(~remove_rows) = zscore(ivT.Relapse(~remove_rows));
        ivZT.Extinction = esc;
        ivZT.Relapse = rel;
        % recall will be nan if the animal did not make a head entry on
        % reinstatement day. these should not be excluded from PCA. 
        rec(~isnan(ivT.Recall) & ~remove_rows) = zscore(ivT.Recall(~isnan(ivT.Recall) & ~remove_rows));
        ivZT.Recall = rec;
    end

    varnames = ivZT.Properties.VariableNames;
    prednames = varnames(varnames ~= "ID" & varnames ~= "Sex" & varnames ~= "Strain" & varnames ~= "Acquire");

    % get Severity scores and classes
    Severity = nansum(ivZT{:, prednames}')';
    Class = cell([height(Severity) 1]);
    Class(Severity>1.5) = {'High'};
    Class(Severity>-1.5 & Severity<1.5) = {'Mid'};
    Class(Severity<-1.5) = {'Low'};
    Class = categorical(Class);
    ivZT.Severity = Severity;
    ivZT.Class = Class;
    
    % remove rows to be excluded from PCA from returned table
    ivZT(remove_rows, :) = [];
    
end
