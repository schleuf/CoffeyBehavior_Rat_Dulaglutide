% Title: main_RatSABehavior_20240429
% Author: Kevin Coffey, Ph.D.
% Affiliation: University of Washington, Psychiatry
% email address: mrcoffey@uw.edu  
% Last revision: 22-May 2024
% LKM revision: 08/14/24, edited rawVariableExtractor to show PR data
% SS revisions: 11/5/24 - 

% ------------- Description --------------
% This is the main analysis script for Golden Oral Fentanyl SA Behavior.
% ----------------------------------------

%% ------------- USER INPUTS --------------
close all
clear all

% IMPORT PATHS
%   main_folder: Path to respository folder
%   beh_datapath: Path to the folder containing all data files to be included in master data table. 
%                 Used to generate a new masterTable if createNewMasterTable == true
%   masterTable_flnm: Path to the masterTable .mat file loaded in if createNewMasterTable == false
%   masterSheet_flnm: Path to key describing information specific to each animal
%   BE_intake_canonical_flnm: Path to sheet recording drug concentration, dose, and intake of 
%                 Behavioral Economics sessions. Only used if runType == 'BE'
%   experimentKey_flnm: % Path to documenting session dates and types for each experiment run
main_folder = pwd;
cd(main_folder)
addpath(genpath(main_folder))
masterTable_flnm = '.\data_masterTable.mat'; 
beh_datapath = {'.\All Behavior'}; 
masterSheet_flnm = '.\Coffey Dula Master Key.xlsx'; 
% BE_intake_canonical_flnm = '.\2024.12.09.BE Intake Canonical.xlsx'; 
experimentKey_flnm = '.\Experiment Key.xlsx'; 

% MISC. SETTINGS
%   runNum: 'all' or desired runs separated by underscores (e.g. '1', '1_3_4', '3_2')
%   runType: 'all', 'ER' (Extinction Reinstatement), 'BE' (Behavioral Economics), 'SA' (Self Administration)
%   createNewMasterTable: If true, generates & saves a new master table from medPC files in datapath. 
%                         If false, reads mT in from masterTable_flnm
%   firstHour: If true, acquire data from the first-hour of data and analyze in addition to the full sessions
%   excludeData: If true, excludes data based on the 'RemoveSession' column of masterSheet
%   acquisition_thresh: To be labeled as "Acquire", animal must achieve an average number of infusions during the 
%                       acquisition_testPeriod greater than this threshold
%   acquisition_testPeriod: Determines sessions to average infusions across before applying acquisition_thresh. 
%                       second value can be 'all', 'first', or 'last'. If 'first' or 'last', there should be a 3rd
%                       value giving the number of days to average across, or it will default to 1.  
%   pAcq:true: plot aquisition histogram to choose threshold

runNum = '2'; 
runType = 'SA'; 
createNewMasterTable = false; 
firstHour = false; 
excludeData = true; 
acquisition_thresh = 10; 
acquisition_testPeriod = {'Training', 'first', 5};
interpWeights = true;
interpWeight_sessions = [0,6,11,16,21];
run_BE_analysis = false;
run_withinSession_analysis = false;
run_individualSusceptibility_analysis = false;

% FIGURE SETTINGS
% Note: if figures are generated they are also saved.
%   saveTabs: If true, save .mat and/or .csv of analyzed datasets
%   pAcq: If true, 
%   dailyFigs: If true, generate daily figures from dailySAFigures.m
%   pubFigs: If true, generate publication figures from pubSAFigures.m
%   indivIntake_figs: If true, generate figures for individual animal behavior across & within sessions
%   groupIntake_figs: If true, generate figures grouped by sex, strain, etc. for animal behavior across & within sessions
%   groupOralFentOutput_figs: If true, generate severity figures
%   figsave_type: Cell of char variables listing all image data types to save figures as
saveTabs = true;
pAcq = true;
dailyFigs = false;
pubFigs = true;
indivIntake_figs = false;
groupIntake_figs = false;
groupOralFentOutput_figs = false;
figsave_type = {'.png', '.pdf'};

% color settings chosen for publication figures. SSnote: haven't been implemented across most figure-generating functions yet. 
gramm_Dula_Sex_colors = {'hue_range',[40 310],'lightness_range',[95 65],'chroma_range',[50 90]};
gramm_Cont_Sex_colors = {'hue_range',[85 -200],'lightness_range',[85 75],'chroma_range',[75 90]};
gramm_Condition_Acq_colors = {'hue_range',[25 385],'lightness_range',[95 60],'chroma_range',[50 70]};
col_M_Dula = [0, 187/255, 144/255];
col_F_Dula = [1, 107/255, 74/255];
col_M_Cont = [163/255, 137/255, 1];
col_F_Cont = [198/255, 151/255, 0];

% SAVE PATHS
% - Each dataset run (determined by runNum and runType) will have its own
%   folder created in the allfig_savefolder. 
% - All other paths will be subfolders of allfig_savefolder designated to 
%   the various figure types and matlab data saved. 
% - Currently only daily & publication figures are saved with current date in
%   the file name, so be aware of overwrite risk for other figures.
allfig_savefolder = 'Outputs\';
dailyfigs_savepath = 'Daily Figures\';
pubfigs_savepath = 'Publication Figures\';
indivIntakefigs_savepath = 'Individual Intake Figures\';
groupIntakefigs_savepath ='Group Intake Figures\'; 
groupOralFentOutput_savepath = 'Severity Output\';
tabs_savepath = 'Behavior Tables\';


%% ------------- HOUSEKEEPING --------------

dt = char(datetime('today')); % Used for Daily & Publication figure savefile names

runNum = categorical(string(runNum));
runType = categorical(string(runType));
if runType == 'all'
    runType = categorical(["ER", "BE", "SA"]);
end

% Import Master Key
opts = detectImportOptions(masterSheet_flnm);
opts = setvartype(opts, {'TagNumber','ID','Cage','Sex','TimeOfBehavior'}, 'categorical'); % Columns of the master key to be pulled in
mKey = readtable(masterSheet_flnm, opts);

% Create subdirectories
toMake = {tabs_savepath, dailyfigs_savepath, pubfigs_savepath, ...
          indivIntakefigs_savepath, groupIntakefigs_savepath, groupOralFentOutput_savepath};
new_dirs = makeSubFolders(allfig_savefolder, runNum, runType, toMake, excludeData, firstHour);
sub_dir = new_dirs{1};
if firstHour
    fH_sub_dir = new_dirs{2};
end

% import experiment key
expKey = readtable(experimentKey_flnm);

%% ------------- IMPORT DATA --------------
if createNewMasterTable
    mT = createMasterTable(beh_datapath, masterSheet_flnm, experimentKey_flnm, 'data_masterTable');
else
    load(masterTable_flnm)
end

%% ------------- FILTER DATA --------------
% exclude data
if excludeData
    mT = removeExcludedData(mT, mKey);
end

% get index for different experiments 
dex = getExperimentIndex(mT, runNum, runType);

% hackymakelifebetterlater
if any(contains(fieldnames(dex), 'ER'))
    if isempty(dex.ER)
        runType(runType == 'ER') = [];
    end
end
if any(contains(fieldnames(dex), 'BE'))
    if isempty(dex.BE)
        runType(runType == 'BE') = [];
    end
end

% Determine Acquire vs Non-acquire
% SS note: this is running on all data in mT, not just the group specified
% for analysis. Plot will reflect ALL data
Acquire = getAcquire(mT, acquisition_thresh, acquisition_testPeriod, pAcq);
if ~any(ismember(mT.Properties.VariableNames, 'Acquire'))
    mT=[mT table(Acquire)];
else
    mT.Acquire = Acquire;
end

% Weight Interpolation
if interpWeights
    mT = interpoweight(mT, interpWeight_sessions);
end

% Get data from the first hour of the session 
if firstHour
    hmT = getFirstHour(mT);
end

%% ------------- GROUP STATISTICS --------------
groupStats = struct;
if firstHour; hour_groupStats = struct; end
for et = 1:length(runType)
    groupStats.(char(runType(et))) = grpstats(mT(dex.(char(runType(et))),:), ["Sex", "Treatment", "Session"], ["mean", "sem"], ...
                          "DataVars",["ActiveLever", "InactiveLever", "EarnedInfusions", "HeadEntries", "Latency", "Intake"]);
    if firstHour
        hour_groupStats.(char(runType(et))) = grpstats(hmT(dex.(char(runType(et))),:),["Sex", "Treatment", "Session"], ["mean", "sem"], ...
                                   "DataVars",["ActiveLever", "InactiveLever", "EarnedInfusions", "HeadEntries", "Latency", "Intake"]);
    end
    if saveTabs
        writeTabs(mT(dex.(char(runType(et))),:), [sub_dir, tabs_savepath, 'run_', char(runNum), '_exp_', char(runType(et)), '_inputData'], {'.mat', '.xlsx'})
        writeTabs(groupStats.(char(runType(et))), [sub_dir, tabs_savepath, 'run_', char(runNum), '_exp_', char(runType(et)), '_GroupStats'], {'.mat', '.xlsx'})
        if firstHour
            writeTabs(hmT(dex.(char(runType(et))),:), [fH_sub_dir, tabs_savepath, 'run_', char(runNum), '_exp_', char(runType(et)), '_inputData'], {'.mat', '.xlsx'})
            writeTabs(hour_groupStats.(char(runType(et))), [fH_sub_dir, tabs_savepath, 'run_', char(runNum), '_exp_', char(runType(et)), '_GroupStats'], {'.mat', '.xlsx'})
        end
    end
end

%% ------------- FIGURES FOR DAILY SPOT CHECKS --------------

if dailyFigs
    %Generate a set of figures to spotcheck data daily
    dailySAFigures(mT, runType, dex, [sub_dir, dailyfigs_savepath], figsave_type);
    % close all
    if firstHour
        dailySAFigures(hmT, runType, dex, [fH_sub_dir, dailyfigs_savepath], figsave_type)
        % close all
    end
end

%% ------------- PUBLICATION-QUALITY FIGURES --------------

if pubFigs %  && strcmp(runType, 'ER')
    pubSAFigures(mT, runType, dex, [sub_dir, pubfigs_savepath], figsave_type);
    if firstHour 
        pubSAFigures(hmT, runType, dex, [fH_sub_dir, pubfigs_savepath], figsave_type); 
    end
    close all;
end

%% ------------- BEHAVIORAL ECONOMICS ANALYSIS --------------
% Note: This section is pulling intake data from '2024.12.09.BE Intake Canonical.xlsx.' 
%       Other analyses and figures in the main analysis script pull intake data from 'Experiment Key.xlsx'
if any(ismember(runType, 'BE')) && run_BE_analysis
    fig_colors = {[.5,.5,.5], col_F_c57, col_M_c57, col_F_CD1, col_M_CD1};
    BE_processes(mT(dex.BE, :), expKey, BE_intake_canonical_flnm, sub_dir, indivIntake_figs, ...
                 groupIntake_figs, saveTabs, fig_colors, indivIntakefigs_savepath, groupIntakefigs_savepath, ...
                 tabs_savepath, figsave_type);
    if firstHour
        BE_processes(hmT(dex.BE, :), expKey, BE_intake_canonical_flnm, fH_sub_dir, indivIntake_figs, ...
                     groupIntake_figs, saveTabs, fig_colors, indivIntakefigs_savepath, groupIntakefigs_savepath, ...
                     tabs_savepath, figsave_type);
    end
end

%% ------------- WITHIN-SESSION ANALYSIS --------------

if run_withinSession_analysis
    fig_colors = {[.5,.5,.5], col_F_c57, col_M_c57, col_F_CD1, col_M_CD1};
    [mTDL, mPressT, mDrugsLT] = WithinSession_Processes(mT, dex, sub_dir, indivIntake_figs, indivIntakefigs_savepath, groupIntake_figs, groupIntakefigs_savepath, saveTabs, tabs_savepath, figsave_type, fig_colors);
end

%% ------------- LINEAR MIXED EFFECTS MODELS --------------

statsname=[sub_dir, tabs_savepath, 'Oral SA Group Stats '];
saveList = {};
% Training
data = mT(mT.sessionType == 'Training',:);
dep_var = ["Intake", "EarnedInfusions", "HeadEntries", "Latency", "ActiveLever", "InactiveLever"];
lme_form = " ~ Sex*Session + (1|TagNumber)";
xlabel('Responses/mg/mL'); % ??? why did kevin put this here
ylabel('Fentanyl Intake (μg/kg)'); % ??? why did kevin put this here

if ~isempty(data)
    Training_LMEstats = getLMEstats(data, dep_var, lme_form);
    if saveTabs
        save([statsname, 'SA'], 'Training_LMEstats');
    end
end

if any(ismember(runType,'ER'))

    % Extinction
    data = mT(mT.sessionType=='Extinction',:);
    dep_var = ["HeadEntries", "Latency", "ActiveLever", "InactiveLever"];
    lme_form = " ~ Sex*Session + (1|TagNumber)";
    if ~isempty(data)
        Extinction_LMEstats = getLMEstats(data, dep_var, lme_form);
    end

    % Reinstatement
    data = mT(mT.sessionType=='Reinstatement',:);
    dep_var = ["HeadEntries", "Latency", "ActiveLever", "InactiveLever"];
    lme_form = " ~ Sex + (1|TagNumber)";
    if ~isempty(data)
        Reinstatement_LMEstats = getLMEstats(data, dep_var, lme_form);
    end
    if saveTabs
        if exist("Extinction_LMEstats", "var")
            save([statsname, 'Extinction'], 'Extinction_LMEstats');
        end
        if exist("Reinstatement_LMEstats", "var")
            save([statsname, 'Reinstatement'], 'Reinstatement_LMEstats');
        end
    end

elseif any(ismember(runType,'BE'))

    % BehavioralEconomics
    data = mT(mT.sessionType=='BehavioralEconomics',:);
    dep_var = ["Intake", "EarnedInfusions", "HeadEntries", "Latency", "ActiveLever", "InactiveLever"];
    lme_form = " ~ Sex + (1|TagNumber)";
    if ~isempty(data)
        BehavioralEconomics_LMEstats = getLMEstats(data, dep_var, lme_form);
        if saveTabs
            save([statsname, 'BE'], 'BehavioralEconomics_LMEstats');
        end
    end
end

%% ------------- INDIVIDUAL SUSCEPTABILITY ANALYSIS --------------
% 1) Calculate individual susceptibility (IS) metrics
% 2) Calculate z-scores for each IS metric, sum these scores for each animal to get the Severity score
% 3) Get correlations between IS metric z-scores & correlation plot, calculated within the following groupings: 
%   - all animals, C57s, CD1s, Males, Females, Male C57s, Female C57s, Male CD1s, Female CD1s
% 4) Make violin plots of IS metrics in the following group pairs:
%   - C57s & CD1s, Males & Females, C57 Males & C57 Females, CD1 Males & CD1 Females
% 5) Generate PCA plots from IS metrics that show all animals against the first 3 and the first 2 principle components. Animals are marked with respect to Strain and Sex.
%       SSnote: add calculation of PCA for subgroups
%
% Individual Susceptibility Metrics
%   1) Intake = total fentanyl consumption in SA (ug/kg)
%   2) Seeking = total head entries in SA
%   3) Cue Association = HE Latency in SA 
%   4) Escalation = Slope of intake in SA
%   5) Extinction = total active lever presses during extinction
%   6) Persistance = slope of extinction active lever presses
%   7) Flexibility = total inactive lever presses during extinction
%   8) Relapse = total presses during reinstatement
%   9) Cue Recall = HE Latency in reinstatement 
% 
% Experiment-dependent use cases: 
%   If the experiment contains ER experiment data, animals involved in these
%       experiments will be analyzed with respect to all IS metrics. 
%   If the experiment contains non-ER experiment data, all animals will be
%       analyzed with respect to Intake, Seeking, Cue Association, and Extinction. 
%
% .mat files saved during the process
%   - ivT: contains IS metrics for all animals
%   - ivZT: contains IS metric Z-scores (separately calculated and saved for ER & nonER groups)
%   - correlations: contains the correlations calculated for all subgroups,
%                   (separately calculated and saved for ER and non ER groups)

if run_individualSusceptibility_analysis

    % subgroups of z-scored data to run correlations across
    corrGroups = {{{'all'}}, ...
                  {{'Strain', 'c57'}}, ...
                  {{'Strain', 'CD1'}}, ...
                  {{'Sex', 'Male'}}, ...
                  {{'Sex', 'Female'}}, ...
                  {{'Strain', 'c57'}, {'Sex', 'Male'}}, ...
                  {{'Strain', 'c57'}, {'Sex', 'Female'}}, ...
                  {{'Strain', 'CD1'}, {'Sex', 'Male'}}, ...
                  {{'Strain', 'CD1'}, {'Sex', 'Female'}}}; 

    % groups to show comparison violin plots for each individual
    % susceptibility metric
    violSubsets = {{'all'}, {'all'}, {'Strain', 'c57'}, {'Strain', 'CD1'}};
    violGroups = {'Strain', 'Sex', 'Sex', 'Sex'};
    violLabels = {'Strain', 'Sex', 'c57 Sex', 'CD1 Sex'};

    pcaGroups = corrGroups;

    ivT = IS_processes(mT, dex, runType, corrGroups, violSubsets, ...
                       violGroups, violLabels, pcaGroups, sub_dir, ...
                       saveTabs, tabs_savepath, groupOralFentOutput_figs, ...
                       groupOralFentOutput_savepath, figsave_type);
    if firstHour
        fH_ivT = IS_processes(hmT, dex, runType, corrGroups, violSubsets, ...
                              violGroups, violLabels, pcaGroups, fH_sub_dir, ...
                              saveTabs, tabs_savepath, groupOralFentOutput_figs, ...
                              groupOralFentOutput_savepath, figsave_type);
    end
end
