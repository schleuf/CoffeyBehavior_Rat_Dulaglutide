function saveFigsByType(f, path, figsave_type)
    for fst = 1:length(figsave_type)
        if strcmp(figsave_type{fst}, '.pdf')
            exportgraphics(f,[path, figsave_type{fst}], 'ContentType','vector')
        else
            exportgraphics(f,[path, figsave_type{fst}]);
        end
        disp(['saved figure: ', path, figsave_type{fst}])
    end
end