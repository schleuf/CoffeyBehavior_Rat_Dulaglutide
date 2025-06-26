function CorrFig(ct, prednames, sub_dir, subfolder, suffix, figsave_type)

    f = figure('Position',[1 1 700 600]);
    imagesc(ct,[-1 1]); % Display correlation matrix as an image
    colormap(flipud(brewermap([],'RdBu')));
    a = colorbar();
    a.Label.String = 'Rho';
    a.Label.FontSize = 12;
    a.FontSize = 12;
    set(gca, 'XTick', 1:length(prednames), 'XTickLabel', prednames, 'XTickLabelRotation', 45, 'FontSize', 12); % set x-axis labels
    set(gca, 'YTick', 1:length(prednames), 'YTickLabel', prednames, 'YTickLabelRotation', 45, 'FontSize', 12); % set x-axis labels
    box off
    set(gca,'LineWidth',1.5,'TickDir','out')
    title(strrep(suffix(2:end), '_', ' '))
    saveFigsByType(f, [sub_dir, subfolder, 'Correlation_table', suffix], figsave_type)
    close(f)
end