function SliceFlipper(optionalVolumePath)
% Interactively flips slices (left-right) from a single multi-page TIFF file.
% Uses keyboard shortcuts for all actions.
% Navigates slices based on the order in _processing_decisions.txt if valid.
% Uses and updates the common _processing_decisions.txt file.
% Can take an optional full path to the TIFF volume as input.

    % --- Initial Setup: File Handling ---
    inputFileFullPath = '';
    if nargin > 0 && ~isempty(optionalVolumePath) && exist(optionalVolumePath, 'file')
        [~, ~, ext] = fileparts(optionalVolumePath);
        if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
            inputFileFullPath = optionalVolumePath;
            disp(['Using provided volume path: ', inputFileFullPath]);
        else
            warning('Provided path is not a TIFF file. Please select a file manually.');
        end
    end

    if isempty(inputFileFullPath)
        [fileName, pathName] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'}, ...
                                         'Select Multi-Page TIFF File');
        if isequal(fileName, 0) || isequal(pathName, 0)
            disp('No file selected. Exiting.');
            return;
        end
        inputFileFullPath = fullfile(pathName, fileName);
    end
    
    try
        tiffInfo = imfinfo(inputFileFullPath);
        numSlices = numel(tiffInfo);
        if numSlices == 0
            errordlg('The selected TIFF file contains no images.', 'File Error');
            return;
        end
    catch ME
        errordlg(['Error reading TIFF file info: ' ME.message], 'File Error');
        return;
    end

    gui_data.originalSliceImages = cell(numSlices, 1);
    hWaitBar = waitbar(0, 'Loading slices from TIFF file...');
    try
        for i = 1:numSlices
            waitbar(i/numSlices, hWaitBar, sprintf('Loading slice %d/%d...', i, numSlices));
            gui_data.originalSliceImages{i} = imread(inputFileFullPath, i, 'Info', tiffInfo);
        end
    catch ME
        if ishandle(hWaitBar); close(hWaitBar); end
        errordlg(['Error loading slices from TIFF: ' ME.message], 'Image Loading Error');
        return;
    end
    if ishandle(hWaitBar); close(hWaitBar); end

    [filePathStr, baseName, ~] = fileparts(inputFileFullPath);
    gui_data.processingDecisionsFilename = fullfile(filePathStr, strcat(baseName, '_processing_decisions.txt'));

    % --- GUI Creation ---
    screenSize = get(0, 'ScreenSize');
    gui_aspect_ratio = 1.6; 
    gui_width_fraction = 0.5;
    gui_width_px = screenSize(3) * gui_width_fraction;
    gui_position = [ ...
        (screenSize(3)-gui_width_px)/2, ... 
        (screenSize(4)-gui_width_px/gui_aspect_ratio)/2, ...
        gui_width_px, gui_width_px/gui_aspect_ratio];

    gui_fig = figure('Name', 'Slice Flipper', ...
        'NumberTitle', 'off', 'Toolbar', 'none', 'Menubar', 'none', ...
        'Color', 'w', 'Units', 'pixels', 'Position', gui_position, ...
        'CloseRequestFcn', @(src,evt) callback_close_gui_request(src), ...
        'KeyPressFcn', @callback_keypress); 

    gui_data.numSlices = numSlices;
    gui_data.currentSequentialPosition = 1; % Position in the display sequence
    gui_data.flipState = zeros(numSlices, 1); % 0 for no flip, 1 for flip L/R
    gui_data.displaySequenceOriginalIndices = (1:numSlices)'; % Default: original order

    % Load existing decisions
    if exist(gui_data.processingDecisionsFilename, 'file') && numSlices > 0
        try
            loadedTable = readtable(gui_data.processingDecisionsFilename, 'Delimiter', '\t', 'ReadVariableNames', true);
            if height(loadedTable) == numSlices && ...
               all(ismember({'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'}, loadedTable.Properties.VariableNames))
                
                % Load FlipState
                flipData = loadedTable.FlipState;
                if isnumeric(flipData)
                    gui_data.flipState = loadedTable.FlipState;
                else 
                    tempFlipState = zeros(numSlices,1);
                    for k_f=1:numSlices
                        if iscell(flipData) && strcmpi(strtrim(string(flipData{k_f})), "1")
                             tempFlipState(k_f) = 1;
                        elseif isstring(flipData) && strcmpi(strtrim(flipData(k_f)), "1")
                             tempFlipState(k_f) = 1;
                        end
                    end
                    gui_data.flipState = tempFlipState;
                    if any(tempFlipState) || ~all(cellfun('isempty', cellstr(string(flipData)))) 
                        warning('FlipState column in %s was not purely numeric. Attempted conversion.', gui_data.processingDecisionsFilename);
                    end
                end
                disp(['Loaded flip states from: ', gui_data.processingDecisionsFilename]);

                % Load and set display sequence based on NewOrderOriginalIndex
                loadedNewOrderCol = loadedTable.NewOrderOriginalIndex;
                isValidPermutation = ~any(isnan(loadedNewOrderCol)) && ...
                                     length(unique(loadedNewOrderCol)) == numSlices && ...
                                     all(ismember(loadedNewOrderCol, 1:numSlices));
                if isValidPermutation
                    gui_data.displaySequenceOriginalIndices = loadedNewOrderCol;
                    disp('Display sequence set from saved order in decisions file.');
                else
                    disp('Saved order in decisions file is not a complete permutation. Using original 1-N display sequence.');
                end
            else
                warning('Dimension mismatch or missing columns in %s. Using defaults.', gui_data.processingDecisionsFilename);
            end
        catch ME_load
            warning('Error loading/parsing decisions from %s: %s. Using defaults.', gui_data.processingDecisionsFilename, ME_load.message);
        end
    end
    
    gui_data.imageAxes = axes('Parent', gui_fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
    gui_data.imageAxes.Colormap = colormap('gray'); 
    axis(gui_data.imageAxes, 'image', 'off'); 
    
    gui_data.imageHandle = imagesc(gui_data.imageAxes, []); 
    gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 10);
    
    guidata(gui_fig, gui_data); 
    
    display_current_slice(gui_fig); 
    
    if gui_data.numSlices > 0 
        uiwait(gui_fig); 
    elseif ishandle(gui_fig) 
        delete(gui_fig);
        disp('No slices to display. GUI closed.');
    end
end

% --- Display Function ---
function display_current_slice(fig)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end

    % Get the original slice index based on the current position in the display sequence
    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentSequentialPosition);
    
    imgData = gui_data.originalSliceImages{originalSliceIdx};

    flip_status_str = 'Normal';
    if gui_data.flipState(originalSliceIdx) == 1 % Flip state is based on original index
        imgData = fliplr(imgData);
        flip_status_str = 'FLIPPED';
    end

    set(gui_data.imageHandle, 'CData', imgData);
    axis(gui_data.imageAxes, 'image'); 
    axis(gui_data.imageAxes, 'off');   

    title_str = {sprintf('Display Pos: %d/%d (Original Slice Index: %d) - %s', ...
                         gui_data.currentSequentialPosition, gui_data.numSlices, originalSliceIdx, flip_status_str), ...
                 'Keys: Left/Right = Change Slice | ''f'' or Ctrl+Arrow = Flip | ''s'' = Save & Close | Esc = Close'};
    set(gui_data.titleHandle, 'String', title_str);
end

% --- Keypress Callback ---
function callback_keypress(fig, eventdata)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end

    ctrl_on = any(strcmp(eventdata.Modifier,'control'));
    originalSliceIdxToFlip = gui_data.displaySequenceOriginalIndices(gui_data.currentSequentialPosition);

    switch eventdata.Key
        case 'leftarrow'
            if ctrl_on
                gui_data.flipState(originalSliceIdxToFlip) = ~gui_data.flipState(originalSliceIdxToFlip);
            else
                gui_data.currentSequentialPosition = max(1, gui_data.currentSequentialPosition - 1);
            end
        case 'rightarrow'
            if ctrl_on
                gui_data.flipState(originalSliceIdxToFlip) = ~gui_data.flipState(originalSliceIdxToFlip);
            else
                gui_data.currentSequentialPosition = min(gui_data.numSlices, gui_data.currentSequentialPosition + 1);
            end
        case 'f'
            gui_data.flipState(originalSliceIdxToFlip) = ~gui_data.flipState(originalSliceIdxToFlip);
        case 's' 
            guidata(fig, gui_data); 
            save_and_close_confirmed(fig);
            return; 
        case 'escape'
            callback_close_gui_request(fig); 
            return; 
    end
    guidata(fig, gui_data);
    display_current_slice(fig);
end

% --- Close GUI Request Callback (for 'x' button or Esc) ---
function callback_close_gui_request(fig)
    gui_data = guidata(fig); 
    if isempty(gui_data) || ~isfield(gui_data, 'numSlices') || gui_data.numSlices == 0
        if ishandle(fig); delete(fig); end 
        return;
    end

    choice = questdlg('Save changes to flip states before closing?', ...
        'Confirm Close', 'Save & Close', 'Discard & Close', 'Cancel', 'Save & Close');
    switch choice
        case 'Save & Close'
            save_and_close_confirmed(fig);
        case 'Discard & Close'
            if ishandle(fig); delete(fig); end
        case 'Cancel'
            return; 
    end
end

% --- Save and Close Action (confirmed) ---
function save_and_close_confirmed(fig)
    gui_data = guidata(fig);
    try
        save_processing_decisions(gui_data);
        disp('Flip decisions saved.');
    catch ME
        errordlg(['Error saving decisions: ', ME.message], 'Save Error');
        return; 
    end
    if ishandle(fig); delete(fig); end
end

% --- Save Processing Decisions Utility ---
function save_processing_decisions(gui_data)
    decisionsMatrix = [];
    if exist(gui_data.processingDecisionsFilename, 'file')
        try
            loadedTable = readtable(gui_data.processingDecisionsFilename, 'Delimiter', '\t', 'ReadVariableNames', true);
            if height(loadedTable) == gui_data.numSlices && ...
               all(ismember({'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'}, loadedTable.Properties.VariableNames))
                
                tempMatrix = zeros(gui_data.numSlices, 3);
                tempMatrix(:,1) = loadedTable.OriginalIndex;
                
                if islogical(loadedTable.FlipState) 
                    tempMatrix(:,2) = double(loadedTable.FlipState);
                elseif isnumeric(loadedTable.FlipState)
                    tempMatrix(:,2) = loadedTable.FlipState;
                else
                     warning('Unexpected data type for FlipState in loaded file. Re-initializing this column.');
                     tempMatrix(:,2) = gui_data.flipState; 
                end

                if isnumeric(loadedTable.NewOrderOriginalIndex)
                    tempMatrix(:,3) = loadedTable.NewOrderOriginalIndex;
                else 
                    tempNewOrder = NaN(gui_data.numSlices,1);
                    if iscell(loadedTable.NewOrderOriginalIndex) || isstring(loadedTable.NewOrderOriginalIndex)
                        for k_no = 1:gui_data.numSlices
                           val_str = string(loadedTable.NewOrderOriginalIndex(k_no)); 
                           if ~strcmpi(strtrim(val_str), "nan") 
                               val = str2double(val_str);
                               if ~isnan(val); tempNewOrder(k_no) = val; end
                           end
                        end
                    end
                    tempMatrix(:,3) = tempNewOrder;
                    if any(~isnan(tempNewOrder)) || ~all(cellfun('isempty', cellstr(string(loadedTable.NewOrderOriginalIndex))))
                        warning('NewOrderOriginalIndex column was not purely numeric. Attempted conversion.');
                    end
                end
                decisionsMatrix = tempMatrix;
            else
                 warning('Dimension or column name mismatch in %s. Initializing.', gui_data.processingDecisionsFilename);
            end
        catch ME_read
            warning('Could not read/parse existing decisions file %s to preserve other columns: %s. Initializing.', ...
                    gui_data.processingDecisionsFilename, ME_read.message);
        end
    end

    if isempty(decisionsMatrix) 
        decisionsMatrix = zeros(gui_data.numSlices, 3);
        if gui_data.numSlices > 0
            decisionsMatrix(:,1) = (1:gui_data.numSlices)';
            decisionsMatrix(:,3) = (1:gui_data.numSlices)'; 
        else 
            decisionsMatrix = zeros(0,3); 
        end
    end
    
    if gui_data.numSlices > 0
        decisionsMatrix(:,2) = gui_data.flipState; 
    end
    
    try
        if gui_data.numSlices > 0
            T = array2table(decisionsMatrix, 'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
            writetable(T, gui_data.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t');
        else 
            T_empty = table('Size', [0,3], 'VariableTypes', {'double', 'double', 'double'}, ...
                            'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
            writetable(T_empty, gui_data.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t'); 
        end
        disp(['Processing decisions (with updated flip states) saved to: ', gui_data.processingDecisionsFilename]);
    catch ME_write
        rethrow(ME_write);
    end
end

