function BE_processes(mT, expKey, BE_intake_canonical_flnm, sub_dir, indivIntake_figs, groupIntake_figs, saveTabs, fig_colors, indivIntakefigs_savepath, groupIntakefigs_savepath, tabs_savepath, figsave_type)
    
    [beT, beiT, BE_IndivFit, BE_GroupFit] = BE_Analysis(mT, expKey, BE_intake_canonical_flnm);
    
    if saveTabs
        writeTabs(beT, [sub_dir, tabs_savepath, 'BE_data'], {'.mat', '.xlsx'})
        writeTabs(beiT, [sub_dir, tabs_savepath, 'BE_allIntake'], {'.mat', '.xlsx'})
        writeTabs(BE_IndivFit, [sub_dir, tabs_savepath, 'BE_indivCurveFit'], {'.mat', '.xlsx'})
        writeTabs(BE_GroupFit.subgroup_curveFits, [sub_dir, tabs_savepath, 'BE_groupCurveFit'], {'.mat', '.xlsx'})
    end
    
    if indivIntake_figs
        BE_CurveFit_Fig(BE_IndivFit, repmat(fig_colors(1),1,height(BE_IndivFit)), [sub_dir, indivIntakefigs_savepath], figsave_type);
    end

    if groupIntake_figs

        % GROUP CURVE FIT FIGS
        BE_CurveFit_Fig(BE_GroupFit.subgroup_curveFits, fig_colors, [sub_dir, groupIntakefigs_savepath], figsave_type);

        % Elasticity and Q0
        figure
        hold on
        for i = 1:height(BE_IndivFit)
            if BE_IndivFit.Sex(i) == 'Female' && BE_IndivFit.Strain(i) == 'c57'
                col = fig_colors{2};
                lab = 'C57 Female';
                a = scatter(BE_IndivFit.Alpha(i), BE_IndivFit.fQ0(i), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
            elseif BE_IndivFit.Sex(i) == 'Male' && BE_IndivFit.Strain(i) == 'c57'
                col = fig_colors{3};
                lab = 'C57 Male';
                b = scatter(BE_IndivFit.Alpha(i), BE_IndivFit.fQ0(i), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
            elseif BE_IndivFit.Sex(i) == 'Female' && BE_IndivFit.Strain(i) == 'CD1'
                col = fig_colors{4};
                lab = 'CD1 Female';
                c = scatter(BE_IndivFit.Alpha(i), BE_IndivFit.fQ0(i), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
            elseif BE_IndivFit.Sex(i) == 'Male' && BE_IndivFit.Strain(i) == 'CD1'
                col = fig_colors{5};
                lab = 'CD1 Male';
                d = scatter(BE_IndivFit.Alpha(i), BE_IndivFit.fQ0(i), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
            end
        end
        l = legend([a, b, c, d], {'C57 Female', 'C57 Male', 'CD1 Female', 'CD1 Male'});
        xlabel('Elasticity')
        ylabel('Q_0')
        hold off


        % GROUP BEHAVIOR FIGS
        xvar = ["unitPrice_log", "unitPrice_log"];
        xlab = ["Fentanyl Concentration (ug/mL)", "Fentanyl Concentration (ug/mL)"];
        xtick = log2([4.5 14.3 45.5]); % SSnote: unhardcode
        xticklab = {'222','70','20'}; % SSnote: unhardcode. also why are we not just plotting against concentration here in the first place
        
        grammOptions = {'color', beT.Sex, 'lightness', beT.Strain, 'subset', 'placeholder'};
        legOptions = {'color', 'Sex', 'lightness', 'Strain'};
        axOptions = {'XTick', xtick, 'XTickLabels', xticklab,'TickDir','Out','LineWidth',1.5};
        pointOptions = {'base_size',9};    
        
        % yvar, ylab, subset, and figname define the different group figures 
        % that will be generated (the rest of the figure parameters are shared
        % across them)
        yvar = {["measuredIntake", "ActiveLever"], ...
                ["Latency", "EarnedInfusions"], ...
                ["measuredIntake", "ActiveLever"], ...
                ["Latency", "EarnedInfusions"]};

        ylab = {["Fentanyl Intake (μg/kg)", "Active Lever Presses"], ...
                ["Head Entry Latency", "Rewards"], ...
                ["Fentanyl Intake (μg/kg)", "Active Lever Presses"], ...
                ["Head Entry Latency", "Rewards"]};
        
        subset = {beT.Acquire=='Acquire' & beT.Strain =='c57', ...
                  beT.Acquire=='Acquire' & beT.Strain =='c57', ...
                  beT.Acquire=='Acquire' & beT.Strain =='CD1', ...
                  beT.Acquire=='Acquire' & beT.Strain =='CD1'};

        figname = {'BE Intake and Active Lever Grouped by Sex for c57', ...
                   'BE Latency and Rewards Grouped by Sex for c57', ...
                   'BE Intake and Active Lever Grouped by Sex for CD1', ...
                   'BE Latency and Rewards Grouped by Sex for CD1'};

        colorOptions = {{'hue_range',[40 310],'lightness_range',[95 65],'chroma_range',[50 90]},...
                        {'hue_range',[40 310],'lightness_range',[95 65],'chroma_range',[50 90]},...
                        {'hue_range',[85 -200],'lightness_range',[85 75],'chroma_range',[75 90]},...
                        {'hue_range',[85 -200],'lightness_range',[85 75],'chroma_range',[75 90]}};
        
        for i = 1:length(yvar)
            grammOptions{end} = subset{i};
            gramm_GroupFig(beT, xvar, yvar{i}, xlab, ylab{i}, [sub_dir, groupIntakefigs_savepath, figname{i}], figsave_type, ...
                           'GrammOptions', grammOptions, 'LegOptions', legOptions, 'PointOptions',pointOptions,'AxOptions', axOptions,'ColorOptions',colorOptions{i});
        end

    end
end