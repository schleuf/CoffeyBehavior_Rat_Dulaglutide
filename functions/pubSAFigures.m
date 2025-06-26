function pubSAFigures(mT, runType, dex, figFold, figsave_type)
% pubSAFigures Generates Clean Subset of Figures for Publication
% pT = the master behavior table from main_MouseSABehavior
% dt = current date time variable
% figFold = Name of the subfolder to save in (relative, folder name is suffiecient

for rt = 1:length(runType)
    expStr = char(string(runType(rt)));

    pT = mT(dex.(string(expStr)),:);
    if runType(rt) == 'ER'
        titles = {'Self-Administration', 'Session 15', 'Extinction', 'Reinstatement'}; % titles of subplots 1-4 for each figure
    elseif runType(rt) == 'BE'
        titles = {'Self-Administration', 'Session 15', 'Behavioral Economics', 'Retraining'};
    elseif runType(rt) == 'SA'
        titles = {'Self-Administration', 'Session 15'};
    end

    %% Active Lever All (grouped by Sex & Acquisition)
    fnum = 1;
    figNames{fnum} = fullfile(figFold,[expStr, '_ActiveLever_SexAcqCollapsed']);
    subset = {{'Acquire', {'Acquire','NonAcquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'ActiveLever';
    yLabs{fnum} = 'Active Lever';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Acquire};
    orderOptions{fnum} = {'color', {'Female', 'Male'}, 'lightness', {'NonAcquire','Acquire'}};
    legendOptions{fnum} = {'x','Sex'};
    colorOptions{fnum} = {'hue_range',[25 385],'lightness_range',[95 60],'chroma_range',[50 70]};
    donut{fnum} = true;
    forceY{fnum} = false;

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});

    %% Active Lever, Acquirers (grouped by Treatment)
    fnum = 2;
    figNames{fnum} = fullfile(figFold,[expStr, '_ActiveLeverAcquire_TreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'ActiveLever';
    yLabs{fnum} = 'Active Lever';
    cVars{fnum} = 'Treatment';
    grammOptions{fnum} = {'lightness',pT.Treatment};
    orderOptions{fnum} = {'lightness',{'Dulaglutide','Vehicle'}};
    legendOptions{fnum} = {'x','Treatment'};
    colorOptions{fnum} = {'hue_range',[40 310],'lightness_range',[85 35],'chroma_range',[30 70]};
    donut{fnum} = false;
    forceY{fnum} = [200 200 200 200];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});
    
    %% Active Lever, Acquirers (grouped by Sex & Treatment)
    fnum = 3;
    figNames{fnum} = fullfile(figFold,[expStr, '_ActiveLeverAcquire_SexTreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'ActiveLever';
    yLabs{fnum} = 'Active Lever';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Treatment};
    orderOptions{fnum} = {'color',{'Female','Male'}, 'lightness', {'Dulaglutide', 'Vehicle'}};
    legendOptions{fnum} = {'x','Sex'};
    colorOptions{fnum} = {'hue_range',[85 -200],'lightness_range',[85 75],'chroma_range',[75 90]};
    donut{fnum} = false;
    forceY{fnum} = [200 200 200 200];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});
    
    

    %% Inactive Lever, Acquirers (grouped by Sex & Treatment)
    fnum = 4;
    figNames{fnum} = fullfile(figFold,[expStr, '_InactiveLeverAcquire_SexTreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'InactiveLever';
    yLabs{fnum} = 'Inactive Lever';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Treatment};
    orderOptions{fnum} = {'color', {'Female','Male'}, 'lightness', {'Dulaglutide', 'Vehicle'}};
    legendOptions{fnum} = {'x','Sex'};
    colorOptions{fnum} = {'hue_range',[40 310],'lightness_range',[85 35],'chroma_range',[30 70]};
    donut{fnum} = false;
    forceY{fnum} = [200 200 200 200];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});
    
    %% Head Entries, Acquirers (grouped by Sex & Treatment)
    fnum = 5;
    figNames{fnum} = fullfile(figFold,[expStr, '_HeadEntriesAcquire_SexTreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'HeadEntries';
    yLabs{fnum} = 'Head Entries';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Treatment};
    orderOptions{fnum} = {'color', {'Female', 'Male'}, 'lightness', {'Dulaglutide', 'Vehicle'}};
    legendOptions{fnum} = {'x', 'Sex'};
    colorOptions{fnum} = {'hue_range', [40 310], 'lightness_range', [95 65], 'chroma_range',[50 90]};
    donut{fnum} = false;
    forceY{fnum} = [600 600 600 600];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});

    %% Latency, Acquirers (grouped by Sex & Treatment)
    fnum = 6;
    figNames{fnum} = fullfile(figFold,[expStr, '_LatencyAcquire_SexTreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'HeadEntries';
    yVars{fnum} = 'Latency';
    yLabs{fnum} = 'Head Entry Latency';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Treatment};
    orderOptions{fnum} = {'color', {'Female','Male'}, 'lightness', {'Dulaglutide', 'Vehicle'}};
    legendOptions{fnum} = {'x', 'Sex'};
    colorOptions{fnum} = {'hue_range', [40 310], 'lightness_range', [95 65], 'chroma_range', [50 90]};
    donut{fnum} = false;
    forceY{fnum} = [200 200 200 200];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});

    %% Intake, Acquirers (grouped by Sex & Treatment)
    fnum = 7;
    figNames{fnum} = fullfile(figFold,[expStr, '_IntakeAcquire_SexTreatCollapsed']);
    subset = {{'Acquire', {'Acquire'}}};
    subInds{fnum} = getSubInds(pT, subset);
    yVars{fnum} = 'Intake';
    yLabs{fnum} = 'Fentanyl Intake (ug/kg)';
    cVars{fnum} = 'Sex';
    grammOptions{fnum} = {'color', pT.Sex, 'lightness', pT.Treatment};
    orderOptions{fnum} = {'color', {'Female','Male'}, 'lightness', {'Dulaglutide', 'Vehicle'}};
    legendOptions{fnum} = {'x','Sex'};
    colorOptions{fnum} = {'hue_range', [40 310], 'lightness_range', [95 65], 'chroma_range', [50 90]};
    donut{fnum} = false;
    forceY{fnum} = [4000 4000 4000 4000];

    plotPubFig(pT, expStr, yVars{fnum}, yLabs{fnum}, forceY{fnum}, cVars{fnum}, subInds{fnum}, titles, figNames{fnum}, figsave_type, donut{fnum}, ...
        'GrammOptions', grammOptions{fnum}, 'ColorOptions', colorOptions{fnum}, ...
        'OrderOptions', orderOptions{fnum}, 'LegendOptions', legendOptions{fnum});

end
end

%%
function [g] = plotPubFig(pT, runType, yVar, yLab, forceY, cVar, subInd, titles, figName, figsave_type, donut, varargin)
    % try
    p = inputParser;
    addParameter(p, 'GrammOptions', {});            
    addParameter(p, 'ColorOptions', {})
    addParameter(p, 'OrderOptions', {});          
    addParameter(p, 'LegendOptions', {});
    
    parse(p, varargin{:});
    
    if strcmp(runType, 'ER')
        sp_subInd = {subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training'), ...
            subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training') & pT.Session>14, ...
            subInd & (pT.sessionType=='Extinction'), ...
            subInd & pT.sessionType=='Reinstatement'};
        xLim = {[0, 15.5], [0.5 2.5], [15.5 25.5], [0.5 2.5]};
    elseif strcmp(runType, 'BE')
        sp_subInd = {subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training'), ...
            subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training') & pT.Session>14, ...
            subInd & (pT.sessionType=='BehavioralEconomics'), ...
            subInd & (pT.sessionType=='ReTraining')};
        xLim = {[0, 15.5], [0.5 2.5], [15.5 20.5], [0.5 2.5]};
    elseif strcmp(runType, 'SA')
        sp_subInd = {subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training'), ...
            subInd & (pT.sessionType=='PreTraining' | pT.sessionType=='Training') & pT.Session>14};
        xLim = {[0, 15.5], [0.5 2.5]};
    end

    stat_set = {{'geom',{'black_errorbar','point','line'},'type','sem','dodge',0,'setylim',1,'width',1}, ...
        {'geom',{'black_errorbar','bar'},'type','sem','dodge',1.75,'width',1.5}};
    point_set = {{'base_size', 9}, {'base_size', 6}};
    if length(titles)==2
        f1=figure('Position',[1 300 600 300]);
    else
        f1=figure('Position',[1 300 1200 300]);
    end
    clear g;

    for sp = 1:length(sp_subInd)
        if any(sp_subInd{sp})
            if mod(sp,2) == 1
                g(1,sp)=gramm('x',pT.Session,'y',pT.(yVar),'subset', sp_subInd{sp}, p.Results.GrammOptions{:});
                g(1,sp).stat_summary(stat_set{1}{:});
                g(1,sp).set_point_options('markers',{'o','s'}, point_set{1}{:});
                g(1,sp).set_names('x','Session','y', yLab,'color',cVar);
                if any(sp_subInd{2})
                    g(1,sp).no_legend;
                end
            elseif mod(sp,2) == 0
                g(1,sp)=gramm('x', pT.(cVar),'y', pT.(yVar), 'subset', sp_subInd{sp}, p.Results.GrammOptions{:});
                g(1,sp).stat_summary(stat_set{2}{:});
                g(1,sp).set_point_options('markers',{'o','s'}, point_set{2}{:});
                g(1,sp).geom_jitter('alpha',.6,'dodge',1.75,'width',0.05);
                g(1,sp).set_names(p.Results.LegendOptions{:});
            end
            g(1,sp).set_text_options('font','Helvetica','base_size',13,'legend_scaling',.75,'legend_title_scaling',.75);
            g(1,sp).set_color_options(p.Results.ColorOptions{:});
            g(1,sp).set_order_options(p.Results.OrderOptions{:});
            g(1,sp).set_title(titles{sp});
        end
    end
    
    g.draw;

    % set limits for the y axis
    if isnumeric(forceY)
       yMax=forceY;
    else
        for sp = 1:length(sp_subInd)
            if any(sp_subInd{sp})
                for s = 1:length(g(1,sp).results.stat_summary)
                    maxStat(s) = nanmax(g(1,sp).results.stat_summary(s).yci(:));
                end
                yMax(sp) = nanmax(maxStat);
            end
        end
        yMax = (ceil(yMax/10)*10)+10;
        if length(yMax) == 2
            yMax(1:2) = max(yMax(1:2));
        end
        if length(yMax) == 4  
            yMax(3:4) = max(yMax(3:4));
        end
    end
    
    for sp = 1:length(sp_subInd)
        if any(sp_subInd{sp})
            %g(1,sp).axe_property('LineWidth', 1.5, 'XLim', xLim{sp}, 'YLim', [0 yMax(sp)], 'TickDir','out');
        
            % Title
            set(g(1,sp).title_axe_handle.Children ,'FontSize',12);
        
            if mod(sp,2) == 1
                % Marker Manipulation
                set(g(1,sp).results.stat_summary(1).point_handle,'MarkerEdgeColor',[0 0 0]);
                set(g(1,sp).results.stat_summary(2).point_handle,'MarkerEdgeColor',[0 0 0]);
                if length(g(1,sp).results.stat_summary) > 2
                    set(g(1,sp).results.stat_summary(3).point_handle,'MarkerEdgeColor',[0 0 0]);
                    set(g(1,sp).results.stat_summary(4).point_handle,'MarkerEdgeColor',[0 0 0]);
                end
            elseif mod(sp,2) == 0
                % Marker Manipulation
                set(g(1,sp).results.geom_jitter_handle(1),'MarkerEdgeColor',[0 0 0]);
                set(g(1,sp).results.geom_jitter_handle(2),'MarkerEdgeColor',[0 0 0]);
                set(g(1,sp).results.stat_summary(1).bar_handle,'EdgeColor',[0 0 0]);
                set(g(1,sp).results.stat_summary(2).bar_handle,'EdgeColor',[0 0 0]);
                if length(g(1,sp).results.stat_summary) > 2
                    set(g(1,sp).results.geom_jitter_handle(3),'MarkerEdgeColor',[0 0 0]);
                    set(g(1,sp).results.geom_jitter_handle(4),'MarkerEdgeColor',[0 0 0]);
                    set(g(1,sp).results.stat_summary(3).bar_handle,'EdgeColor',[0 0 0]);
                    set(g(1,sp).results.stat_summary(4).bar_handle,'EdgeColor',[0 0 0])
                end
            end
        end
    end

    % Axes Limits
    set(g(1,1).facet_axes_handles,'YLim',[0 yMax(1)],'XLim',[0 15.5]);
    if any(sp_subInd{2})
        if strcmp(cVar,'Treatment')
            set(g(1,2).facet_axes_handles,'YLim',[0 yMax(2)]);
        elseif strcmp(cVar,'Sex')
            set(g(1,2).facet_axes_handles,'YLim',[0 yMax(2)],'XTickLabel',{char(9792),char(9794)});
        end
    end

     % Axes Limits
    if length(sp_subInd) > 2
        if strcmp(runType, 'ER') && any(sp_subInd{3})
            set(g(1,3).facet_axes_handles,'YLim',[0 yMax(3)],'XLim',[15 25.5],'XTick',[16 20 25],'XTickLabel',{'1','5','10'});
        elseif strcmp(runType, 'BE') && any(sp_subInd{3})
            set(g(1,3).facet_axes_handles,'YLim',[0 yMax(3)],'XLim',[15 20.5],'XTick',[16 18 20],'XTickLabel',{'222','70','20'});
        end
        
        if strcmp(cVar,'Treatment') && any(sp_subInd{4})
            set(g(1,4).facet_axes_handles,'YLim',[0 yMax(4)]);
        elseif strcmp(cVar,'Sex')  && any(sp_subInd{4})
            set(g(1,4).facet_axes_handles,'YLim',[0 yMax(4)],'XTickLabel',{char(9792),char(9794)});
        end
    end
    
    % Remove & Move Axes
    if any(sp_subInd{2})
        set(g(1,2).facet_axes_handles,'YColor',[1 1 1]);
        set(g(1,2).facet_axes_handles,'YLabel',[],'YTick',[]);
    
        if height(g(1, 2).results.geom_jitter_handle)==4
            pos1=g(1,2).facet_axes_handles.OuterPosition;
            set(g(1,2).facet_axes_handles,'OuterPosition',[pos1(1)-.04,pos1(2),pos1(3)-.04,pos1(4)]);
            pos2=g(1,2).title_axe_handle.OuterPosition;
            set(g(1,2).title_axe_handle,'OuterPosition',[pos2(1)-.05,pos2(2),pos2(3),pos2(4)]);
        elseif height(g(1, 2).results.geom_jitter_handle)==2
            pos1=g(1,2).facet_axes_handles.OuterPosition;
            set(g(1,2).facet_axes_handles,'OuterPosition',[pos1(1)-.03,pos1(2),pos1(3)-.1,pos1(4)]);
            pos2=g(1,2).title_axe_handle.OuterPosition;
            set(g(1,2).title_axe_handle,'OuterPosition',[pos2(1)-.08,pos2(2),pos2(3),pos2(4)])
        end
    end

    if length(sp_subInd) > 2
        % Remove & Move Axes
        if any(sp_subInd{3})
            pos5=g(1,3).facet_axes_handles.OuterPosition;
            set(g(1,3).facet_axes_handles,'OuterPosition',[pos5(1),pos5(2),pos5(3)-.065,pos5(4)]);
            pos6=g(1,3).title_axe_handle.OuterPosition;
            set(g(1,3).title_axe_handle,'OuterPosition',[pos6(1)-.035,pos6(2),pos6(3),pos6(4)]);
        end
        if any(sp_subInd{4})
            set(g(1,4).facet_axes_handles,'YColor',[1 1 1]);
            set(g(1,4).facet_axes_handles,'YLabel',[],'YTick',[]);
            if height(g(1, 2).results.geom_jitter_handle)==4
                pos3=g(1,4).facet_axes_handles.OuterPosition;
                set(g(1,4).facet_axes_handles,'OuterPosition',[pos3(1)-.09,pos3(2),pos3(3)-.04,pos3(4)]);
                pos4=g(1,4).title_axe_handle.OuterPosition;
                set(g(1,4).title_axe_handle,'OuterPosition',[pos4(1)-.09,pos4(2),pos4(3),pos4(4)]);
            elseif height(g(1, 2).results.geom_jitter_handle)==2
                pos3=g(1,4).facet_axes_handles.OuterPosition;
                set(g(1,4).facet_axes_handles,'OuterPosition',[pos3(1)-.08,pos3(2),pos3(3)-.11,pos3(4)]);
                pos4=g(1,4).title_axe_handle.OuterPosition;
                set(g(1,4).title_axe_handle,'OuterPosition',[pos4(1)-.135,pos4(2),pos4(3),pos4(4)]);
            end
        end
    end

     for sp = 1:length(sp_subInd)
        if any(sp_subInd{sp})
            set(g(1,sp).facet_axes_handles,'LineWidth', 1.5, 'TickDir','out');
        end
     end
    
    % Export Figure
    for fst = 1:length(figsave_type)
        if strcmp(figsave_type{fst}, '.pdf')
            exportgraphics(gcf,[figName, figsave_type{fst}], 'ContentType','vector')
        elseif strcmp(figsave_type{fst}, '.fig')
            savefig(gcf,[figName, figsave_type{fst}]);
        else
            exportgraphics(gcf,[figName, figsave_type{fst}]);
        end
        disp(['saved figure: ', figName, figsave_type{fst}])
    end
    
    
    if donut
        plotDonut(pT, subInd, g, figName, figsave_type)
    end
    % catch
    % disp('oh no! had to skip a figure, wonder why??')
    % end
end
    
function plotDonut(pT, subInd, g, figName, figsave_type)
    % Donut Chart for Overlay
    groupStats = grpstats(pT(pT.Session==1 & subInd, :),["Sex","Treatment","Acquire"],["mean","sem"],"DataVars",'ActiveLever');
    % groupStats = sortrows(groupStats,"Acquire",'descend');
    f2=figure('Position',[1 300 575 575]);
    d=donutchart(groupStats.GroupCount, strrep(groupStats.Properties.RowNames, '_', ' '));
    d.InnerRadius=.45;
    % SS note: can't figure out how to make sure these colors correspond
    % to the figure when there are groups with no data
    % for s = 1:length(g(1,1).results.stat_summary)
    %     d.ColorOrder(s,1:3)=g(1,1).results.stat_summary(s).point_handle.MarkerFaceColor; % SS note: the hell is going on here w/ the stat summary indexing?
    % end
    d.FontSize=12.5;
    d.FontName='Arial Rounded MT Bold';
    
    % Export Figure
    for fst = 1:length(figsave_type)
        if strcmp(figsave_type{fst}, '.pdf')
            exportgraphics(gcf,[figName, '_Donut', figsave_type{fst}], 'ContentType','vector')
        elseif strcmp(figsave_type{fst}, '.fig')
            savefig(gcf,[figName, figsave_type{fst}]);
        else
            exportgraphics(gcf,[figName, '_Donut', figsave_type{fst}]);
        end
        disp(['saved figure: ', figName, '_Donut', figsave_type{fst}])
    end

end