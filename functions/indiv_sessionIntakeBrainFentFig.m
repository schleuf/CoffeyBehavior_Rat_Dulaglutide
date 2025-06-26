function indiv_sessionIntakeBrainFentFig(xvar, yvar, figpath, figsave_type)
    f=figure('Position',[100 100 400 800]);
    clear g
    g(1,1)=gramm('x',xvar{1}, 'y', yvar{1}); 
    g(1,1).stat_bin('normalization','cumcount','geom','stairs','edges',0:1:180);
    g(2,1)=gramm('x',xvar{2}, 'y', yvar{2}); % SS note: what's this *1000 for?
    g(2,1).geom_line();
    g(1,1).axe_property('LineWidth',1.5,'FontSize',13,'XLim',[0 180],'tickdir','out');
    g(1,1).set_names('x','Session Time (m)','y','Cumulative Infusions');
    g(2,1).axe_property('LineWidth',1.5,'FontSize',13,'XLim',[0 180],'tickdir','out');
    g(2,1).set_names('x','Session Time (m)','y','Brain Fentanyl Concentration ug/kg');
    g.draw();
    saveFigsByType(f, figpath, figsave_type)
    close(f);
end