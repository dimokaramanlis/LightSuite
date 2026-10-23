function settings = parseSettingsFile(filePath)
%parseSettingsFile Parses a text configuration file into a MATLAB struct.
%
%   SETTINGS = parseSettingsFile(FILEPATH) reads the configuration file
%   specified by FILEPATH. Each line in the file containing an '=' is
%   treated as a key-value pair.
%
%   The function returns a struct, SETTINGS, where each field corresponds
%   to a key from the file.
%
%   - Lines starting with '%' or '#' are treated as comments and ignored.
%   - Empty lines are ignored.
%   - If a setting is not specified in the file, a default value is used.
%   - If the file does not exist, it is created with all default values so
%     the user can edit it, and the defaults are returned.
%   - Values are evaluated as MATLAB expressions, allowing for numbers,
%     strings (in single quotes), logicals (true/false), and arrays.
%
%   Example file format:
%       % Main settings
%       slicethickness = 80
%       regchan = 'cy3'
%       use_gpu = true
%       roi_coords = [10 10 200 200]
%

    % --- 1. Define Default Settings ---
    % These values will be used if they are not specified in the text file.
    % This makes your analysis scripts robust to missing settings.
    disp('Loading default settings...');
    settings.denoisedapi = false;
    settings.slicethickness = 150;
    settings.px_process = 5.0;
    settings.px_register = 20.0;
    settings.px_atlas = 10.0;
    settings.atlasaplims = [180 1079];
    settings.regchan = 'dapi';
    settings.medianfiltreg = false;
    settings.use_gpu = false;


    % --- 2. Check if the settings file exists ---
    % If not, write the defaults there so the user has a file to edit.
    if ~exist(filePath, 'file')
        try
            writeSettingsFile(filePath, settings);
            warning(['Settings file not found. Created %s with default values.\n' ...
                'Edit it to match your data and re-run.'], filePath);
        catch ME
            warning('Settings file not found: %s and could not be created (%s). Using default values.', ...
                filePath, ME.message);
        end
        return; % Exit the function, returning the default struct
    end
    
    disp(['Found settings file. Reading overrides from: ' filePath]);

    % --- 3. Open and read the file ---
    % Use a try-catch block to ensure the file is always closed properly,
    % even if an error occurs during parsing.
    fileID = fopen(filePath, 'r');
    if fileID == -1
        warning('Could not open settings file: %s. Using all default values.', filePath);
        return;
    end
    
    % Ensure file gets closed when the function exits or errors
    cleanupObj = onCleanup(@() fclose(fileID));

    % --- 4. Parse the file line by line ---
    while ~feof(fileID)
        line = strtrim(fgetl(fileID)); % Read a line and remove whitespace

        % Skip empty lines or lines that are comments
        if isempty(line) || startsWith(line, '%') || startsWith(line, '#')
            continue;
        end

        % Check for a key-value pair (must contain '=')
        if contains(line, '=')
            parts = strsplit(line, '=', 'CollapseDelimiters', false);
            if numel(parts) >= 2
                key = strtrim(parts{1});
                % Join the rest of the parts in case the value contains '='
                valueStr = strtrim(strjoin(parts(2:end), '=')); 

                % Ensure the key is a valid MATLAB struct field name
                key = matlab.lang.makeValidName(key);

                try
                    % IMPORTANT: eval() is used to convert the string value
                    % into a MATLAB type (e.g., number, logical, string, array).
                    % This is very powerful but assumes the settings file is
                    % from a trusted source and not malicious.
                    evaluatedValue = eval(valueStr);
                    
                    % Overwrite the default setting with the value from the file
                    settings.(key) = evaluatedValue;
                    fprintf('  - Overriding setting: %s = %s\n', key, valueStr);
                    
                catch ME
                    warning('Could not parse value for key "%s". Error: %s. Skipping.', key, ME.message);
                end
            end
        end
    end
    
    disp('Finished parsing settings file.');

end

function writeSettingsFile(filePath, settings)
% Write SETTINGS as 'key = value' lines that parseSettingsFile can read back.
    fileID = fopen(filePath, 'w');
    if fileID == -1
        error('Could not open %s for writing.', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fileID));

    fprintf(fileID, '%% Local Settings File - %s\n', char(datetime('now')));
    fprintf(fileID, '%% Edit the values below. Lines starting with %% or # are comments.\n');
    fprintf(fileID, '%% Settings not found in this file will use internal defaults.\n\n');
    for fn = fieldnames(settings)'
        val = settings.(fn{1});
        if ischar(val) || isstring(val)
            valstr = ['''' strrep(char(val), '''', '''''') ''''];
        else
            valstr = mat2str(val);
        end
        fprintf(fileID, '%s = %s\n', fn{1}, valstr);
    end
end
