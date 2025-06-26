function [f] = plotViolins(ivT, yVars, yLabs, group, colorlabel)
    clear g
    f = figure('position',[100 100 width(yVars)*180 300]);
    numDat = length(ivT.Intake); 
    x = nan([1,numDat]);
    groupsets = unique(ivT.(group));
    
    if length(groupsets) > 1
        x(ivT.(group) == categorical(groupsets(1))) = 1;
        x(ivT.(group) == categorical(groupsets(2))) = 1; 
    
        for y = 1:length(yVars)
            g(1,y)=gramm('x',x,'y',ivT.(yVars{y}),'color',ivT.(group));
            if strcmp('Strain',colorlabel)
            g(1,y).set_color_options('hue_range',[25 385],'lightness_range',[95 60],'chroma_range',[50 70],'n_color',2);
            g(1,y).set_order_options('color', groupsets)
            elseif strcmp('Sex',colorlabel)
            g(1,y).set_color_options('hue_range',[25 385],'lightness_range',[95 60],'chroma_range',[50 70],'n_color',2);
            g(1,y).set_order_options('color', groupsets)
            elseif strcmp('c57 Sex',colorlabel)
            g(1,y).set_color_options('hue_range',[40 310],'lightness_range',[95 65],'chroma_range',[50 90],'n_color',2);
            g(1,y).set_order_options('color',{'Female','Male'})
            elseif strcmp('CD1 Sex',colorlabel)
            g(1,y).set_color_options('hue_range',[85 -200],'lightness_range',[85 75],'chroma_range',[75 90],'n_color',2);
            g(1,y).set_order_options('color',{'Female','Male'})
            end
            g(1,y).stat_violin('normalization', 'width', 'fill', 'transparent','extra_y', 0, 'half', 1); %'extra_y', 0, 'half', 1, 
            g(1,y).geom_jitter('width',.05,'dodge',.5,'alpha',.75);
            g(1,y).axe_property('LineWidth',1.5,'FontSize',12,'Font','Helvetica','XLim',[0.5 1.5],'TickDir','out'); %'YLim',[0 1200]
            g(1,y).set_names('x','','y', yLabs{y},'color', '');
            g(1,y).set_point_options('base_size', 6);
            if y~=length(yVars)
            g(1,y).no_legend();
            end
            % g(1,y).set_title('');
        end
        
         g.draw;
    
        for i=1:width(g)
           g(1,i).results.geom_jitter_handle(1).MarkerEdgeColor = [0 0 0]; 
           g(1,i).results.geom_jitter_handle(2).MarkerEdgeColor = [0 0 0];
        end

    else
        disp(['only one value for grouping "', group, '", skipping violin plots...'])
    end

end