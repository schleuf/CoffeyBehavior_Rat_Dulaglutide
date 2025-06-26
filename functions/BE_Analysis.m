function [beT, beiT, BE_IndivFit, BE_GroupFit] = BE_Analysis(mT, expKey, BE_intake_canonical_flnm)

    beT=mT(mT.sessionType == 'BehavioralEconomics' & mT.Acquire == 'Acquire', :); % Initialize Master Table 
    beT.ID = beT.TagNumber;
    ID=unique(beT.TagNumber);  
    BE_sess = unique(expKey.Session(strcmp(expKey.SessionType,'BehavioralEconomics')));

    % Import Dose and Measured Intake Data
    opts = detectImportOptions(BE_intake_canonical_flnm); 
    beiT=readtable(BE_intake_canonical_flnm, opts);
    beiT.TagNumber = categorical(beiT.TagNumber);
    beiT = beiT(ismember(beiT.TagNumber, beT.TagNumber),:);

    % remove sessions with missing data from beiT
    missing_days = cell([length(ID), 1]);
    for i = 1:length(ID)
        this_sess = beT.Session(beT.TagNumber == ID(i));
        md = find(~ismember(BE_sess, this_sess));
        missing_days{i} = md;
        remove_ind = logical((beiT.TagNumber==ID(i)) .* ismember(beiT.Day, missing_days{i}));
        beiT(remove_ind,:) = [];
    end
    unitPrice=(1000./beiT.Dose__g_ml_);
    unitPrice_log = log2(unitPrice);
    beT = [beT, table(beiT.Day, beiT.measuredIntake, beiT.Dose__g_ml_, unitPrice, unitPrice_log)];
    beT = renamevars(beT, ["Var1", "Var2", "Var3", "Var4", "Var5"], ["Day", "measuredIntake", "Dose", "unitPrice", "unitPrice_log"]);

    % Curve Fit Each Animals Intake over Dose
    [Sex, Strain] = deal(categorical(nan([length(ID), 1])));
    [Alpha, Beta, Elastic, fQ0, knee_x] = deal(nan([length(ID), 1]));
    [fitY, fitX, modX, modY] = deal(cell([length(ID), 1]));

    myfittype = fittype('log(b)*(exp(1)^(-1*a*x))', ...
                        'dependent', {'y'}, 'independent', {'x'}, ...
                        'coefficients', {'a','b'});
    coef_start = [.003, 20000]; 
    coef_lower = [.0003 100];
    coef_upper = [.03 100000];

    for i=1:height(ID)
        Sex(i,1)=unique(beT.Sex(beT.TagNumber==ID(i)));
        Strain(i,1)=unique(beT.Strain(beT.TagNumber==ID(i)));
    
        in=beT.measuredIntake(beT.TagNumber==ID(i));
        price=beT.unitPrice(beT.TagNumber==ID(i));
        price(in == 0) = [];
        in(in == 0) = []; 
        % price=price-min(price)+1; Don't think this is necessary...

        fitX{i} = price;
        fitY{i} = log(in);
        
        if height(in) > 1
            f=fit(fitX{i}, fitY{i}, myfittype,'StartPoint', coef_start, 'lower', coef_lower, 'upper', coef_upper);
            Alpha(i)=f.a;
            Beta(i) = f.b; 
            fQ0(i) = exp(f(1)); 
            modX{i} = 1:50;
            modY{i} = f(modX{i});
            [res_x, idx_x]=knee_pt(log(1:50),f(1:50));
            knee_x(i) = res_x;
            Elastic(i) = idx_x;
        end  
    end

    BE_IndivFit = table(ID, Sex, Strain, Alpha, Beta, Elastic, fQ0, knee_x, fitX, fitY, modX, modY);

    % SSnote: this is terribly redundant but i don't care right now
    subgroups = {ones([height(beT), 1]); ...
                 (beT.Strain == 'c57' & beT.Sex == 'Female'); ...
                 (beT.Strain == 'c57' & beT.Sex == 'Male'); ...
                 (beT.Strain == 'CD1' & beT.Sex == 'Female'); ...
                 (beT.Strain == 'CD1' & beT.Sex == 'Male')};

    ID = ["All"; "C57_Female"; "C57_Male"; "CD1_Female"; "CD1_Male"];

    BE_GroupFit = struct;
    [fitX, fitY, semY, modX, modY] =  deal(mat2cell(nan([length(subgroups), 1]), ones([1,length(subgroups)]), 1));
    [Alpha, Beta, Elastic, fQ0, knee_x] = deal(nan([length(subgroups), 1]));
    
    BE_GroupFit.LME = fitlme(beT,'measuredIntake ~ Concentration + (1|ID)');
    BE_GroupFit.ANOVA = anova(BE_GroupFit.LME,'DFMethod','satterthwaite'); % SSnote: not doing anything with these anovas, should look at them somewhere...
    
    for sg = 1:length(subgroups)      
        intake = beT.measuredIntake(find(subgroups{sg}));
        price = beT.unitPrice(find(subgroups{sg}));
        price(intake == 0) = [];
        intake(intake == 0) = []; 
        price = price-min(price)+1;
        if length(intake) > 1

            uni_price = unique(price);
            mean_intake = arrayfun(@(x) mean(log(intake(price==x))), uni_price);
            sem_intake = arrayfun(@(x) std(log(intake(price==x)))/sqrt(length(intake(price==x))), uni_price);
            fitX{sg} = uni_price;
            fitY{sg} = mean_intake;
            semY{sg} = sem_intake;
    
            % Group BE Demand Curve
            ff=fit(fitX{sg}, fitY{sg}, myfittype, 'StartPoint', coef_start,'lower', coef_lower,'upper',coef_upper);
    
            Alpha(sg) =ff.a;
            Beta(sg) = ff.b;
            fQ0(sg) = exp(ff(1)); 
            modX{sg} = 1:50;
            modY{sg} = ff(modX{sg});
            [res_x, idx_x] = knee_pt(log(1:500),ff(1:500)); % SSnote: why log here, is that interacting weird with plotting later
            knee_x(sg) = res_x;
            Elastic(sg) = idx_x; % SSnote: come back to this
        end
    end

    BE_GroupFit.subgroup_curveFits = table(ID, subgroups, Alpha, Beta, Elastic, fitX, fitY, semY, modX, modY, knee_x);
end