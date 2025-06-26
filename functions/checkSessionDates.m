%% checkSessionDates - 
% SS 2025
function [mT] = checkSessionDates(mT, mKey, expKey, correctFiles, savename)

    %% EXPERIMENT TYPE DEFINITIONS FOR LOGICAL INDEXING
    % need to updating the indexing from this to exclude things that include
    % these session types but also include unstated session types. works for now. 
    exp_types = dictionary(["ER", "BE"], ...
                            {{'SelfAdministration', 'Extinction', 'Reinstatement'}, ...
                             {'SelfAdministration', 'BehavioralEconomics'}});
    
    
    %% PULL SPECIFIC RUN & EXPERIMENT TYPE TO CHECK

    mKey = setExperiment_mKey(mKey, exp_types);

    mT = setExperiment_mT(mT, mKey);

    date = mT.Date;
    sess = mT.Session;
    run = mT.Run;

    if any(ismember(mT.Properties.VariableNames, 'sessionType'))
        type = mT.sessionType;
    else
        type = repmat(categorical("undefined"), [length(sess),1]);
    end
    
    % avoid bugs if session type is undefined
    type(isundefined(type)) = categorical("undefined");
    
    
    %% FIND FILES THAT NEED TO BE CORRECTED
    incorrect_session = [];
    incorrect_type = [];
    
    for d = 1:length(expKey.Date)
        key_date = expKey.Date(d);
        key_session = expKey.Session(d);
        key_type = expKey.SessionType(d);
        key_exp = expKey.Experiment(d);
        key_run = expKey.Run(d);

        wrong_session = find(date== key_date & key_run == run & sess ~= key_session & strcmp(mT.Experiment, char(key_exp)));
        wrong_type = find(date == key_date & key_run == run & type ~= key_type & strcmp(mT.Experiment, char(key_exp)));
    
        incorrect_session = [incorrect_session; wrong_session];
        incorrect_type = [incorrect_type; wrong_type];
        
        if ~isempty(wrong_session)
            disp(['Incorrect Sessions logged for ', char(key_date), '(Session = ', num2str(key_session), ') :'])
            for is = 1:length(wrong_session)
                disp(['     Tag ', char(mT.TagNumber(wrong_session(is))), ' labeled as Session ', num2str(sess(wrong_session(is)))]);
            end
            disp(' ')
        end
    
        if ~isempty(wrong_type)
            disp(['Incorrect Session Types logged for ', char(key_date), '(Session Type = ', key_type{1}, ') :'])
            for is = 1:length(wrong_type)
                disp(['     Tag ', char(mT.TagNumber(wrong_type(is))), ' labeled as Session Type ', char(type(wrong_type(is)))]);
            end
            disp(' ')
        end
    
    end
    
    %% PLOT DATE AGAINST SESSION # BY SESSION TYPE
    figure
    hold on
    
    uniType = unique(type);
    leg = uniType;
    
    for uni = 1:length(uniType)
        ind = find(type == uniType(uni));
        scatter(date(ind), sess(ind));
    end
    
    if ~isempty(incorrect_session)
        scatter(date(incorrect_session), sess(incorrect_session), '+k')
        leg = [leg; categorical("incorrect session #")];
    end
    if ~isempty(incorrect_type)
        scatter(date(incorrect_type), sess(incorrect_type), 'xk')
        leg = [leg; categorical("incorrect session type")];
    end
    
    ylim([nanmin(sess) - 1, nanmax(sess) + ((nanmax(sess)-nanmin(sess))/2)])
    legend(leg, 'Location', 'northwest')
    xlabel('Date')
    ylabel('Session #')
    
    hold off
    
    %% CORRECT MEDPC FILES
    if correctFiles
        if ~isempty(incorrect_session) || ~isempty(incorrect_type)
            mT = corrections(mT, mKey, incorrect_session, incorrect_type, expKey);
        else
            disp('No corrections to be made!')
            disp(' ')
        end
       % save new table
       main_folder = pwd;
       save([main_folder, '\', savename], 'mT')
       disp('saved updated masterTable with Experiment column')
       disp(' ')
    end
    
end


%%
function [mT] = corrections(mT, mKey, incorrect_session, incorrect_type, sdKey)
    
    if ~isempty(incorrect_session)
        % correct incorrect sessions in medPC data files
        mT = correctSessions(mT, mKey, incorrect_session, sdKey);
    end
    
   if ~isempty(incorrect_type)
       % correct incorrect session types in master table
       mT = correctTypes(mT, mKey, incorrect_type, sdKey);
   end
end


function [mT] = correctSessions(mT, mKey, incorrect_session, sdKey)
    strStart = 15;
    strEnd = 20;
    dates = datetime(sdKey.Date);
    template = '    15:       ';
    disp('updating incorrect session numbers in medPC files...')
    for is = 1:length(incorrect_session)
        this_date = mT.Date(incorrect_session(is));
        this_file = mT.FileName{incorrect_session(is)};
        this_tag = mT.TagNumber(incorrect_session(is));
        mKey_ind = (mKey.TagNumber==this_tag);
        this_exp = mKey.Experiment(mKey_ind);
        this_exp = this_exp{1};
        this_run = mKey.Run(mKey_ind);
        

        sdKey_ind = find((dates == this_date) .* (sdKey.Run == this_run) .* strcmp(sdKey.Experiment, this_exp));
        corr_sess = sdKey.Session(sdKey_ind);
        char_sess = [num2str(corr_sess), '.000'];
        if corr_sess/10 < 1
            char_sess = [' ', char_sess];
        end
        fileContent = fileread(this_file);
        lines = splitlines(fileContent);
        
        % can't hardcode the line to change bc different number of blank lines in some files
        chng_line = 0;
        ln = 1;
        while chng_line == 0
            if contains(lines{ln}, template)
                chng_line = ln;
            else
                ln = ln + 1;
            end
        end

        tmp = lines{chng_line};
        tmp(strStart:strEnd) = char_sess;
        lines{chng_line} = tmp;
        
        fid = fopen(this_file, 'w');
        for i = 1:length(lines)
            fprintf(fid, '%s\n', lines{i});
        end
        fclose(fid);
        disp(['updated file: ' this_file])
        
        % update slideSession variable
        sessionType = mT.sessionType(incorrect_session(is));
        if sessionType == 'PreTraining'
            mT.slideSession(incorrect_session(is)) = mT.Session(incorrect_session(is));
        elseif sessionType == 'Training'
            mT.slideSession(incorrect_session(is)) = mT.Session(incorrect_session(is)) + 1;
        elseif sessionType == 'Extinction' || sessionType == 'BehavioralEconomics'
            mT.slideSession(incorrect_session(is)) = mT.Session(incorrect_session(is)) + 2;
        elseif sessionType == 'Reinstatement' || sessionType == 'ReTraining'
            mT.slideSession(incorrect_session(is)) = mT.Session(incorrect_session(is)) + 3;
        else
            mT.slideSession(incorrect_session(is)) = mT.Session(incorrect_session(is));
        end
    end    
    disp(' ')


end

function [mT] = correctTypes(mT, mKey, incorrect_type, sdKey)
    disp('updating incorrect session types in masterTable...')
    for it = 1:length(incorrect_type)
        this_date = mT.Date(incorrect_type(it));
        this_tag = mT.TagNumber(incorrect_type(it));
        mKey_ind = (mKey.TagNumber==this_tag);
        this_exp = mKey.Experiment{mKey_ind};
        this_run = mKey.Run(mKey_ind);

        sdKey_ind = find(strcmp(sdKey.Date,string(this_date)) .* (sdKey.Run == this_run) .* strcmp(sdKey.Experiment, this_exp));
        corr_type = categorical(string(sdKey.SessionType{sdKey_ind}));
        mT.sessionType(incorrect_type(it)) = corr_type;
    end
    disp(['updated session types for ', num2str(length(incorrect_type)), ' sessions'])
    disp(' ')
end


function [cat1_cat2s] = checkUnique(cat1, var1, cat2, var2)
    uni_cat1 = unique(var1);
    cat1_cat2s = cell([length(uni_cat1), 1]);
    for uni = 1:length(uni_cat1)
        this_uni_cat2 = unique(var2(var1 == uni_cat1(uni)));
        cat1_cat2s{uni} = this_uni_cat2;
        num_uni_cat2 = height(this_uni_cat2);

        if num_uni_cat2 > 1
            disp(strcat(string(num2str(num_uni_cat2)), " ", string(cat2), "s for ", string(cat1), " ", string(uni_cat1(uni)), ":"))
            for ud = 1:num_uni_cat2
                disp(['     ', char(categorical(this_uni_cat2(ud)))])
            end
            disp('  ')
            disp('  ')
        end
    end

end


function [mKey] = setExperiment_mKey(mKey, exp_types)
    
    Experiment = cell([length(mKey.Run), 1]);
    dict_keys = keys(exp_types);
    for exp = 1:length(dict_keys)
        sessTypes = exp_types(dict_keys(exp));
        sessTypes = sessTypes{1};
        inds = ones([length(mKey.Run), 1]);

        for st = 1:length(sessTypes)
            inds = inds .* mKey.(string(sessTypes{st}));
        end

        Experiment(find(inds)) = {dict_keys(exp)};
    end
    mKey = [mKey, table(Experiment)];
end


function [mT] = setExperiment_mT(mT, mKey)
    Experiment = cell([length(mKey.Run), 1]);
    for tn = 1:length(mKey.TagNumber)
        this_exp = mKey.Experiment(tn);
        this_exp = this_exp{1};
        tn_inds = (mT.TagNumber == mKey.TagNumber(tn));
        Experiment(find(tn_inds)) = {char(this_exp)};
    end
    if any(ismember(mT.Properties.VariableNames, 'Experiment'))
        mT = removevars(mT,{'Experiment'});
    end
    mT = [mT, table(Experiment)];
end