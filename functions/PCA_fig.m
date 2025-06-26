function PCA_fig(ivZT, PCA, sub_dir, subfolder, suffix, figsave_type)
    
    male_c57 = [0, 187/255, 144/255];
    female_c57 = [1, 107/255, 74/255];
    male_CD1 = [163/255, 137/255, 1];
    female_CD1 = [198/255, 151/255, 0];

    groups = fieldnames(PCA);

    for p = 1:length(groups)
        coeff = PCA.(groups{p}).coeff;
        score = PCA.(groups{p}).score;
        prednames = PCA.(groups{p}).prednames;
        ivZT_inds = PCA.(groups{p}).ivZT_inds;
        this_ivZT = ivZT(ivZT_inds,:);


        % 3D figure
        f1=figure('color', 'w', 'position', [100 100 800 650]);
        h1 = biplot(coeff(:,1:3), 'Scores', score(:,1:3),...
            'Color', 'b', 'Marker', 'o', 'VarLabels', prednames);

        % set metric vectors' appearance
        for i = 1:length(prednames) 
            h1(i).Color=[.5 .5 .5];    
            h1(i).LineWidth=1.5;
            h1(i).LineStyle=':';
            h1(i).MarkerSize=4;
            h1(i).MarkerFaceColor=[.0 .0 .0];
            h1(i).MarkerEdgeColor=[0 .0 0];
            h1(i).Marker='o';
        end
    
        % remove extra line objects (not sure why these exist) 
        for i = length(prednames) + 1 : length(prednames) * 2
            h1(i).Marker='none';
        end
    
        % format text for metric vector labels
        for i = 1 + (length(prednames) * 2) : 3 * length(prednames)
            h1(i).FontSize = 11;
            h1(i).FontWeight = 'bold';
        end
        data_ind1 = length(h1) - height(this_ivZT);
        R = rescale(this_ivZT.Severity,4,18);
        for i=data_ind1:length(h1)-1

            h1(i).MarkerEdgeColor=[0 .0 0];
            h1(i).MarkerSize=R(i-data_ind1 + 1);
            if this_ivZT.Strain(i - data_ind1 + 1) == 'c57' && this_ivZT.Sex(i - data_ind1 + 1) == 'Male'
                h1(i).MarkerFaceColor = male_c57;
            elseif this_ivZT.Strain(i - data_ind1 + 1) == 'c57' && this_ivZT.Sex(i - data_ind1 + 1) == 'Female'
                h1(i).MarkerFaceColor = female_c57;
            elseif this_ivZT.Strain(i - data_ind1 + 1) == 'CD1' && this_ivZT.Sex(i - data_ind1 + 1) == 'Male'
                 h1(i).MarkerFaceColor = male_CD1;
            elseif this_ivZT.Strain(i - data_ind1 + 1) == 'CD1' && this_ivZT.Sex(i - data_ind1 + 1) == 'Female'
                 h1(i).MarkerFaceColor = female_CD1;
            end
            
    
        end
        pbaspect([1,1,1]);
        set(gca,'LineWidth',1.5,'TickDir','in','FontSize',14);
        grid off
        title(strrep(groups{p}, '_', ' '));

        saveas(f1, [sub_dir, subfolder '3D_PCA_Vectors_', groups{p}]);
        

        % 2D Figure
        PC1=score(:,1);
        PC2=score(:,2);    
    
        pcTable = [this_ivZT, table(PC1, PC2)];
        f1 = figure('color','w','position',[100 100 800 650]);
        g = gramm('x', pcTable.PC1, 'y', pcTable.PC2, 'color', pcTable.Sex, 'marker', pcTable.Strain, 'lightness', pcTable.Class);
        g.set_order_options('lightness', {'Low', 'Mid', 'High'});
        g.set_color_options('hue_range',[50 542.5],'chroma',80,'lightness',60,'n_color',2);
        g.geom_point();
        g.set_names('x','PC1','y','PC2','color','Sex', 'marker', 'Strain', 'lightness', 'Class');
        g.axe_property('FontSize',12,'LineWidth',1.5,'TickDir','out');
        g.set_order_options('lightness',{'High','Mid','Low'});
        g.set_point_options('base_size',8);
        g.draw;
        for i = 1:height(g.results.geom_point_handle)
           g.results.geom_point_handle(i);
           g.results.geom_point_handle(i).MarkerEdgeColor = [0 0 0];
           g.results.size(i);
        end
        title(strrep(groups{p}, '_', ' '));
        saveFigsByType(f1, [sub_dir, subfolder, '2D_PCA_', groups{p}], figsave_type);
        close(f1)
    end
end