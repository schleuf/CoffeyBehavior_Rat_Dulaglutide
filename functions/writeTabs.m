function writeTabs(data, flnm, fltypes)
    for ft = 1:length(fltypes)
        if strcmp(fltypes(ft), '.mat')
            % MatLab .mat
            save(flnm, 'data');
            varNm = data.Properties.VariableNames;
        elseif strcmp(fltypes(ft), '.xlsx')
            % Excel .xlsx (remove any cell data from table)
            remove_inds = arrayfun(@(x) (class(data.(varNm{x}))=="cell"), 1:length(varNm));
            sub_data = removevars(data,varNm(remove_inds));
            writetable(sub_data, [flnm, '.xlsx'], 'Sheet', 1);
        else
            disp(['Table format "', ftips(ft), '" not recognized, table not saved.'])
        end
    end
end