function [g] = gramm_GroupSubplot(xvar, yvar, xlab, ylab, leg, varargin)

    p = inputParser;
    addParameter(p, 'GrammOptions', {});            
    addParameter(p, 'OrderOptions', {});           
    addParameter(p, 'LegOptions', {});
    addParameter(p, 'ColorOptions', {});
    addParameter(p, 'StatOptions', {});
    addParameter(p, 'AxOptions', {});
    addParameter(p, 'PointOptions', {});
    addParameter(p, 'WrapOptions', {});

    parse(p, varargin{1}{:});
    
    % SSnote: can move these default setting checks to a subfunction
    if isempty(p.Results.ColorOptions)
        ColorOptions = {'hue_range',[0 360],'lightness_range',[85 35],'chroma_range',[30 70]};
    else
        ColorOptions = p.Results.ColorOptions;
    end

    if isempty(p.Results.StatOptions)
        StatOptions = {'geom',{'black_errorbar','point','line'},'type','sem','dodge',0,'setylim',1,'width',1};
    else
        StatOptions = p.Results.StatOptions;
    end

    if isempty(p.Results.AxOptions)
        AxOptions = {'LineWidth', 1.5, 'Font', 10, 'tickdir','out'};
    else
        AxOptions = p.Results.AxOptions;
    end

    g(1,1)=gramm('x',xvar,'y',yvar, p.Results.GrammOptions{:});
    if ~isempty(p.Results.WrapOptions)
        g.facet_wrap(p.Results.WrapOptions{:});
    end
    g(1,1).set_color_options(ColorOptions{:});
    g(1,1).set_point_options(p.Results.PointOptions{:});
    g(1,1).set_text_options('font','Helvetica','base_size',13,'legend_scaling',.75,'legend_title_scaling',.75);
    g(1,1).axe_property(AxOptions{:});
    g(1,1).set_names('x',xlab,'y',ylab, p.Results.LegOptions{:});  

    str_inds = arrayfun(@(x) ischar(StatOptions{x}), 1:length(StatOptions));
    if any(ismember(StatOptions(str_inds), 'cumcount'))
        g(1,1).stat_bin(StatOptions{:});
    else
        g(1,1).stat_summary(StatOptions{:});
    end
    
    if ~leg
        g(1,1).no_legend();
    end
end