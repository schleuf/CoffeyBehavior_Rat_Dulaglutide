function indiv_allSessionFig(tab, subset, xvar, xlab, yvar, ylab, figtitle, facetwrap, figpath, figsave_type, stat_type)
    x = tab.(xvar);
    y = tab.(yvar);
    fw = tab.(facetwrap);
    if strcmp(xvar,'adj_rewLP')
        x = x/60;
    end

    f=figure('Position',[1 1 1920 1080]);
    g=gramm('x', x(subset), 'y', y(subset));
    g.set_color_options('hue_range',[-65 265],'chroma',80,'lightness',70,'n_color',2);
    g.facet_wrap(fw(subset),'scale','independent','ncols',5,'force_ticks',1);
    if strcmp(stat_type, 'cumbin')
        g.stat_bin('normalization','cumcount','geom','stairs','edges',0:1:180);
    else
        g.stat_summary('geom', stat_type,'setylim',1);
    end
    g.axe_property('LineWidth',1.5,'FontSize',12,'XLim',[0 180],'tickdir','out');
    g.set_names('x', xlab,'y',ylab);
    g.set_title(figtitle);
    g.draw;
    for i=1:length(g.facet_axes_handles)
        g.facet_axes_handles(i).Title.FontSize=12;
        set(g.facet_axes_handles(i),'XTick',[0 90 180]);
    end
    saveFigsByType(f, figpath, figsave_type)
    close(f)
end