function tab = loadTableHDF5(loadpath)
    info = h5info(loadpath); % Get structure of HDF5 file
    num_groups = height(info.Groups); % Get number of datasets
    tab = table(); % Initialize empty table

    % Loop through each dataset in the file
    disp(['loading in data from ', loadpath, '...'])
    for i = 1:num_groups
        % get data, column name, and data type to output
        varname = info.Groups(i).Name(2:end); % Get dataset name (column name)
        data = h5read(loadpath, ['/' varname '/data']); % Read dataset
        og_class = h5read(loadpath, ['/' varname '/class']);
        
        % Restore the correct data type
        if width(data) == 1
            if strcmp(og_class, "double")
                data = double(data);
            elseif strcmp(og_class, "categorical")
                data = categorical(data);
            elseif strcmp(og_class, "string")
                data = string(data);
            elseif strcmp(og_class, "char")
                cell_data = cell([height(data), 1]);
                for row = 1:height(data)
                    cell_data{row} = char(data(row));
                end
                data = cell_data;
            elseif strcmp(og_class, "datetime")
                data = datetime(data);
            else
                warning(['loadTableHDF5 not equipped to handle data of type "' char(class(data)), '", skipping column "', char(varname), '...']);
                continue
            end
        else
            cell_data = cell([height(data), 1]);
            if strcmp(og_class, "double")
                for row = 1:height(data)
                
                    last_ind = find(~isnan(data(row,:)), 1, 'last');
                    if ~isempty(last_ind)
                        cell_data{row} = data(row,1:last_ind);
                    else
                        cell_data{row} = [];
                    end
                    cell_data{row} = data(row,1:last_ind)';
                end
               
                data = cell_data;
            else
                warning(['loadTableHDF5 not equipped to handle array data of type "' char(class(data)), '", skipping column "', char(varname), '...'])
                continue
                
            end
        end

        % Add restored data to the table
        tab.(varname) = data;
    end
    disp('HDF5 file loaded successfully');
end