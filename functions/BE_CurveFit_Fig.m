function BE_CurveFit_Fig(fit_tab, colors, figpath, figsave_type)
    yMax = max(arrayfun(@(x) max(fit_tab.modY{x}), 1:height(fit_tab))) + .05;
    yMin = min(arrayfun(@(x) min(fit_tab.modY{x}), 1:height(fit_tab))) - .05;
    fit_lw = 1.5; % line width (not applied to error bars or marker edge)
    font_size = 14; % font size (everything but legend)
    leg_text_scale = .75; % scale of legend font relative to the rest
    leg_line_length = 10;
    marker_size = 75;
    marker_lw = 1;
    error_lw = marker_lw;
    error_capsize = 3;
    ax_lw = 1.5;

    for i = 1:height(fit_tab)
        if ~isempty(fit_tab.modY{i})

            f=figure('Position',[100 100 285 350],'Color',[1 1 1]);
            hold on;
            % plot(log(fit_tab.modX{i}), exp(fit_tab.modY{i}), 'LineWidth', fit_lw, 'Color', colors{i});
            plot(log(fit_tab.modX{i}), fit_tab.modY{i}, 'LineWidth', fit_lw, 'Color', colors{i});
            % e = errorbar(log(fit_tab.fitX{i}), fit_tab.fitY{i}, fit_tab.semY{i}, 'k', 'LineStyle', "none");
            try
            e = errorbar(log(fit_tab.fitX{i}), exp(fit_tab.fitY{i}), exp(fit_tab.semY{i}), 'k', 'LineStyle', "none");
            e.LineWidth = error_lw;
            e.CapSize = error_capsize;
            catch
                disp('Individual, No Error to Make Bars')
            end
            p = scatter(log(fit_tab.fitX{i}), fit_tab.fitY{i}, marker_size, colors{i}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', marker_lw);
            % p = scatter(log(fit_tab.fitX{i}), exp(fit_tab.fitY{i}), marker_size, colors{i}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', marker_lw);
            % plot([log(fit_tab.knee_x(i)), log(fit_tab.knee_x(i))], [min(fit_tab.modY{i}) max(fit_tab.modY{i})],'--k');
            % xlim([0 4.25]);
            ylim([yMin, yMax])
            % ylim([0, 4000])
            xlim([1, log(50)])
            set(gca,'LineWidth',1.5,'tickdir','out','box',0);
            xt=log([4 5 10 25 50]);
            yt=round(exp([6.5 7 7.5 8]));
            set(gca,'XTick',[xt],'XTickLabels',{'4' '5' '10' '25' '50'});
            set(gca,'YTick',[6.5 7 7.5 8],'YTickLabels',yt);
            xlabel('Responses/mg/mL');
            ylabel('Fentanyl Intake (μg/kg)');
            % title(char(strrep(fit_tab.ID(i), '_', ' ')));
            
            % Update all text elements in the figure
            set(findall(f, '-property', 'FontSize'), 'FontSize', font_size);
            set(findall(f, '-property', 'FontName'), 'FontName', 'Helvetica');
            lgd = legend(p, [' ', char(strrep(fit_tab.ID(i), '_', ' '))]);
            try
            lgd.IconColumnWidth = leg_line_length;
            lgd.FontSize = font_size * leg_text_scale;
            catch
            end

            legend('boxoff')
            saveFigsByType(f, [figpath, 'BE_curvefit_', char(fit_tab.ID(i))], figsave_type);

            close(f)
        else
            disp(['no fit data for TagNumber ' char(fit_tab.ID(i)), 'skipping curve-fit plot...']);
        end
    end
end