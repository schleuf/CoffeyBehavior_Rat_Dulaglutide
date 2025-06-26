function ViolinFig(ivT, group, label, includeER, sub_dir, groupOralFentOutput_savepath, figsave_type, violLabels)
    % Violin plots
    

    if includeER
        yVars = {'Intake', 'Seeking', 'Association', 'Escalation', 'Extinction', 'Relapse', 'Recall', 'Severity'};
        yLabs = {' Fentanyl Intake (mg/kg)', 'Seeking (Head Entries)', 'Association (Latency)', ...
                 'Escalation (slope Training Intake)', 'Extinction Responses', 'Relapse (Reinstatement Responses)', 'Recall (Reinstatement Latency)', 'Severity'};
    else
        yVars = {'Intake', 'Seeking', 'Association', 'Escalation',  'Severity'};
        yLabs = {' Fentanyl Intake (mg/kg)', 'Seeking (Head Entries)', 'Association (Latency)', ...
                 'Escalation (slope Training Intake)', 'Severity' };
    end

    f = plotViolins(ivT, yVars, yLabs, group, violLabels);
    saveFigsByType(f, [sub_dir, groupOralFentOutput_savepath, 'Violin', label], figsave_type)
    close(f)
end