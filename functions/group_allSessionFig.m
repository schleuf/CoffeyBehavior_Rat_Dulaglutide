% NOT CURRENTLY USING, KEEPING TEMPORARILY FOR REFERENCE
function group_allSessionFig(tab, subset, xvar, xlab, yvar, ylab, colorGroup, lightGroup, facetwrap, figtitle, figpath, stat_type, figsave_type)
    x = tab.(xvar);
    y = tab.(yvar); 
    cg = tab.(colorGroup);
    lg = tab.(lightGroup);
    if strcmp(xvar,'adj_rewLP')
        x = x/60;
    end
    
    f = figure('units','normalized','outerposition',[0 0 .5 1]);
    g=gramm('x',x(subset),'y',y(subset),'color',cg(subset),'lightness',lg(subset));
    g.set_color_options('hue_range',[50 542.5],'chroma',80,'lightness',60,'n_color',2);
    if ~strcmp(facetwrap,'none')
        fw = tab.(facetwrap);
        g.facet_wrap(fw(subset),'scale','independent','ncols',3,'force_ticks',1,'column_labels',1);
        g.set_names('column','Session');
    end
    
    if strcmp(stat_type, 'cumbin')
        g.stat_bin('normalization','cumcount','geom','stairs','edges',0:1:180);
    else
        g.stat_summary('geom', stat_type,'setylim',1);
    end
    
    g.axe_property('LineWidth',1.5,'FontSize',10,'XLim',[0 180],'tickdir','out');
    g.set_names('x', xlab, 'y',ylab, 'color', colorGroup, 'lightness', lightGroup);
    g.set_title(figtitle);
    g.draw;
    
    yMax = 0;
    
    if strcmp(stat_type, 'cumbin')
        for ss = 1:length(g.results.stat_bin)
            yMax = max([yMax;g.results.stat_bin(ss).counts]);
        end
    else
        for ss = 1:length(g.results.stat_summary)
            if strcmp(stat_type, 'area')
                yMax = max([yMax;g.results.stat_summary(ss).yci(:,2)]);
            elseif strcmp(stat_type, 'line')
                yMax = max([yMax;g.results.stat_summary(ss).y]);
            end
        end
    end
    yMax = yMax + (.05 * yMax);

    for i=1:length(g.facet_axes_handles)
        g.facet_axes_handles(i).YLim=[0 yMax];
    end

    saveFigsByType(f, figpath, figsave_type)
    close(f)
end