function [LME_stats] = getLMEstats(data, dep_var, lme_form)
    LME_stats = struct;
    for dv = 1:length(dep_var)
        LME_stats.(strcat(dep_var(dv), "LME")) = fitlme(data, strcat(dep_var(dv), lme_form));
        LME_stats.(strcat(dep_var(dv), "F")) = anova(LME_stats.(strcat(dep_var(dv), "LME")) ,'DFMethod','satterthwaite');
    end
end