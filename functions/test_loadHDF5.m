
masterTable_flnm = '.\data_masterTable.mat'; % the masterTable .mat file loaded in if createNewMasterTable == false

load(masterTable_flnm)

test_path = '.\masterDataset.hdf5';
saveTableHDF5(mT, test_path);
disp(' ')
newb = loadTableHDF5(test_path);
disp(' ')

%%
colnames = mT.Properties.VariableNames;
loadcolnames = newb.Properties.VariableNames;
problems = false;
if length(colnames) ~= length(loadcolnames)
    warning('whyyyy')
    problems = true;
end
for cn = 1:length(colnames)
    d1 = mT.(colnames{cn});
    d2 = newb.(colnames{cn});
    if ~strcmp(class(d1), class(d2))
        disp(colnames{cn})
        disp('whyyyy')
        problems = true;
    end
    if ~iscell(d1)
        if ~all(d1 == d2)
            z = find(d1~=d2);
            if ~(all(isnan(d1(z))) && all(isnan(d2(z))))
                disp(colnames{cn})
                disp('whyyyy')
                problems = true;
            end
        end
    else
        for row = 1:height(mT)
            if ~all(d1{row} == d2{row})
                disp(colnames{cn})
                disp('whyyyy')
                problems = true;
            end
        end
    end
end

if ~problems
    disp('saved and loaded data are a perfect match!')
else
    disp('you got problems')
end