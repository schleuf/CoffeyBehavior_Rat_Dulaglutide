function saveTableHDF5(tab, savepath)

    % Delete existing file to avoid conflicts
    if isfile(savepath)
        delete(savepath);
    end
    
    % Create an empty HDF5 file before the loop
    fid = H5F.create(savepath);
    H5F.close(fid);
    
    % Loop through each column, creating groups and writing datasets
    disp(['saving data to ', savepath, '...'])
    for i = 1:width(tab)
        % Get column name and data
        varname = tab.Properties.VariableNames{i}; % Get column name
        data = tab{:, i}; % Extract column data
        data_path = ['/' varname '/data']; % Create group path
        class_path = ['/' varname '/class'];
        data = tab.(varname);
        if ~iscell(data)
            og_class = string(class(data));
        else
            og_class = string(class(data{1}));
        end
        % Get datatype for the hdf5 dataset, update the current data type
        % and/or shape as necessary
        if isnumeric(data) % Numeric column
            datatype = "double";
        elseif iscategorical(data) || isstring(data) || isdatetime(data)
            data = string(data);
            datatype = "string";
        elseif iscell(data) && all(cellfun(@isnumeric, data)) % Column with numeric arrays
            max_len = max(cellfun(@numel, data)); % Find max array length
            padded_data = NaN(height(tab), max_len); % Pad shorter arrays with NaNs
            for j = 1:height(tab)
                len = numel(data{j});
                padded_data(j, 1:len) = data{j}; % Fill matrix
            end
            data = padded_data;
            datatype = "double";
        elseif iscell(data) && all(cellfun(@ischar, data)) % Column with numeric arrays
            str_array = strings(size(data));
            for j = 1:height(tab)
                str_array(j) = string(data{j}); % Fill matrix
            end
            data = str_array;
            datatype = "string";
        else
            warning(['saveTableHDF5 is not equipped to handle datatype of column "', varname, '", skipping...'])
            continue
        end
    
        % Create the group for this column
        h5create(savepath, data_path, size(data), 'Datatype', datatype, 'Deflate', 9, 'Chunksize', [100, width(data)]);
        h5create(savepath, class_path, [1,1], 'Datatype', "string");
        h5write(savepath, data_path, data);
        h5write(savepath, class_path, og_class);
    end
    disp('HDF5 file saved successfully');
end