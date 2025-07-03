function SliceOrderEditor(optionalVolumePath)
% Interactively reorders slices using a direct "Move to Position" workflow.
% This is the final, redesigned version for clarity and ease of use.

    % --- Initial Setup & GUI Creation (No changes here) ---
    inputFileFullPath = '';
    if nargin > 0 && ~isempty(optionalVolumePath) && exist(optionalVolumePath, 'file')
        [~, ~, ext] = fileparts(optionalVolumePath);
        if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
            inputFileFullPath = optionalVolumePath;
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
        if numSlices == 0, errordlg('The selected TIFF file contains no images.', 'File Error'); return; end
    catch ME
        errordlg(['Error reading TIFF file info: ' ME.message], 'File Error');
        return;
    end
    gui_data.originalSliceImages = cell(numSlices, 1);
    hWaitBar = waitbar(0, 'Loading slices...');
    try
        for i = 1:numSlices
            waitbar(i/numSlices, hWaitBar);
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
    screenSize = get(0, 'ScreenSize');
    gui_aspect_ratio = 1.6; gui_width_fraction = 0.5;
    gui_width_px = screenSize(3) * gui_width_fraction;
    gui_position = [(screenSize(3)-gui_width_px)/2, (screenSize(4)-gui_width_px/gui_aspect_ratio)/2, gui_width_px, gui_width_px/gui_aspect_ratio];
    gui_fig = figure('Name', 'Slice Order Editor', 'NumberTitle', 'off', 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', 'Units', 'pixels', 'Position', gui_position, 'CloseRequestFcn', @(src,evt) callback_close_gui_request(src), 'KeyPressFcn', @callback_keypress);
    gui_data.numSlices = numSlices;
    gui_data.currentDisplayPosition = 1; 
    gui_data.flipState = zeros(numSlices, 1); 
    gui_data.displaySequenceOriginalIndices = (1:numSlices)'; 
    load_processing_decisions(); 
    gui_data.imageAxes = axes('Parent', gui_fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
    gui_data.imageAxes.Colormap = colormap('gray'); axis(gui_data.imageAxes, 'image', 'off'); 
    gui_data.imageHandle = image(gui_data.imageAxes, []); 
    gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 10);
    gui_data.excludeMarkerHandle = [];
    gui_data.orderTextHandle = text(gui_data.imageAxes, 0, 0, '', 'FontSize', 24, 'Color', 'yellow', 'FontWeight', 'bold', 'BackgroundColor', [0 0 0 0.5], 'VerticalAlignment', 'top');
    guidata(gui_fig, gui_data); 
    display_current_slice(gui_fig); 
    if gui_data.numSlices > 0, uiwait(gui_fig); 
    elseif ishandle(gui_fig), delete(gui_fig); disp('No slices to display. GUI closed.'); end

    function load_processing_decisions()
        if exist(gui_data.processingDecisionsFilename, 'file') && gui_data.numSlices > 0
            try
                loadedTable = readtable(gui_data.processingDecisionsFilename, 'Delimiter', '\t', 'ReadVariableNames', true);
                if height(loadedTable) ~= gui_data.numSlices || ~all(ismember({'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'}, loadedTable.Properties.VariableNames))
                    warning('Decision file has mismatching dimensions or columns. Using defaults.'); return;
                end
                flipData = loadedTable.FlipState;
                if isnumeric(flipData)
                    gui_data.flipState = flipData;
                else 
                    tempFlipState = zeros(gui_data.numSlices,1);
                    for k_f = 1:gui_data.numSlices
                        val_str = strtrim(string(flipData{k_f}));
                        if strcmpi(val_str, "1"), tempFlipState(k_f) = 1;
                        elseif strcmpi(val_str, "-1"), tempFlipState(k_f) = -1;
                        end
                    end
                    gui_data.flipState = tempFlipState;
                end
                loadedOrder = loadedTable.NewOrderOriginalIndex;
                isValidPermutation = ~any(isnan(loadedOrder)) && numel(unique(loadedOrder)) == gui_data.numSlices;
                if isValidPermutation
                    gui_data.displaySequenceOriginalIndices = loadedOrder;
                    disp('Loaded saved order from decisions file.');
                else
                    disp('Saved order is incomplete or invalid. Using default 1-N order.');
                end
            catch ME_load
                warning('Error loading/parsing decisions file: %s. Using defaults.', ME_load.message);
            end
        end
    end
end

function display_current_slice(fig)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    if ishandle(gui_data.excludeMarkerHandle); delete(gui_data.excludeMarkerHandle); gui_data.excludeMarkerHandle = []; end
    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    imgData = gui_data.originalSliceImages{originalSliceIdx};
    currentFlipState = gui_data.flipState(originalSliceIdx);
    status_str = 'Normal';
    if currentFlipState == 1, imgData = fliplr(imgData); status_str = 'FLIPPED';
    elseif currentFlipState == -1, status_str = 'EXCLUDED'; end
    set(gui_data.imageHandle, 'CData', imgData); axis(gui_data.imageAxes, 'image', 'off');   
    if currentFlipState == -1
        [h, w, ~] = size(imgData);
        gui_data.excludeMarkerHandle = text(gui_data.imageAxes, w/2, h/2, 'X', 'FontSize', h/4, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
    end
    orderNumStr = sprintf('%d', gui_data.currentDisplayPosition);
    set(gui_data.orderTextHandle, 'Position', [0.02*size(imgData,2), 0.02*size(imgData,1)], 'String', orderNumStr);
    title_str = {sprintf('Slice at Order Position: %d/%d (Original Index: %d) - %s', ...
                         gui_data.currentDisplayPosition, gui_data.numSlices, originalSliceIdx, status_str), ...
                 'Keys: Nav = Left/Right | Reorder = Enter | Save = s'};
    set(gui_data.titleHandle, 'String', title_str);
    guidata(fig, gui_data);
end

function callback_keypress(fig, eventdata)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    switch eventdata.Key
        case 'leftarrow'
            gui_data.currentDisplayPosition = max(1, gui_data.currentDisplayPosition - 1);
        case 'rightarrow'
            gui_data.currentDisplayPosition = min(gui_data.numSlices, gui_data.currentDisplayPosition + 1);
        case 'f'
            if gui_data.flipState(originalSliceIdx) ~= -1, gui_data.flipState(originalSliceIdx) = ~gui_data.flipState(originalSliceIdx); end
        case 'o'
            if gui_data.flipState(originalSliceIdx) == -1, gui_data.flipState(originalSliceIdx) = 0;
            else, gui_data.flipState(originalSliceIdx) = -1; end
        case {'return', 'enter'}
            reorder_slice_callback(fig); return; 
        case 's', guidata(fig, gui_data); save_processing_decisions(gui_data);return; 
        case 'escape', callback_close_gui_request(fig); return; 
    end
    guidata(fig, gui_data);
    display_current_slice(fig);
end

function reorder_slice_callback(fig)
    gui_data = guidata(fig);
    current_pos = gui_data.currentDisplayPosition;
    slice_to_move_idx = gui_data.displaySequenceOriginalIndices(current_pos);
    
    % --- REDESIGNED WORKFLOW ---
    prompt = {sprintf('Move slice (Orig. Idx: %d) to new order position (1-%d):', slice_to_move_idx, gui_data.numSlices)};
    dlg_title = 'Move to Position';
    answer = inputdlg(prompt, dlg_title, [1 60], {num2str(current_pos)});
    
    if isempty(answer), return; end 
    
    new_pos = round(str2double(answer{1}));
    
    if isnan(new_pos) || new_pos < 1 || new_pos > gui_data.numSlices
        warndlg('Invalid input. Please enter a valid position number.', 'Input Error');
        return;
    end
    
    % --- Call the new, simpler, direct logic ---
    new_order_vec = perform_move_to_position(gui_data.displaySequenceOriginalIndices, slice_to_move_idx, new_pos);
    
    gui_data.displaySequenceOriginalIndices = new_order_vec;
    gui_data.currentDisplayPosition = new_pos; % The new position is exactly what the user entered
    
    guidata(fig, gui_data);
    display_current_slice(fig);
end

function new_order_vec = perform_move_to_position(order_vec, slice_to_move_idx, new_pos)
    % This is the new, simpler logic: remove a slice and insert it at a specific index.
    
    % Remove the slice from the list
    temp_order_vec = order_vec(order_vec ~= slice_to_move_idx);
    
    if new_pos == 1
        % Insert at the beginning
        new_order_vec = [slice_to_move_idx; temp_order_vec];
    elseif new_pos == numel(order_vec)
        % Insert at the end
        new_order_vec = [temp_order_vec; slice_to_move_idx];
    else
        % Insert in the middle
        new_order_vec = [temp_order_vec(1:new_pos-1); ...
                         slice_to_move_idx; ...
                         temp_order_vec(new_pos:end)];
    end
    disp(sprintf('Moved slice %d to position %d.', slice_to_move_idx, new_pos));
end

function callback_close_gui_request(fig)
    choice = questdlg('Save changes to order and states before closing?', 'Confirm Close', 'Save & Close', 'Discard & Close', 'Cancel', 'Save & Close');
    if strcmp(choice, 'Save & Close'), save_and_close_confirmed(fig);
    elseif strcmp(choice, 'Discard & Close'), if ishandle(fig); delete(fig); end
    end
end

function save_and_close_confirmed(fig)
    gui_data = guidata(fig);
    try
        save_processing_decisions(gui_data);
        disp('Order and state decisions saved.');
    catch ME
        errordlg(['Error saving decisions: ', ME.message], 'Save Error');
        return; 
    end
    if ishandle(fig); delete(fig); end
end

function save_processing_decisions(gui_data)
    decisionsMatrix = zeros(gui_data.numSlices, 3);
    decisionsMatrix(:,1) = (1:gui_data.numSlices)'; 
    decisionsMatrix(:,2) = gui_data.flipState;
    decisionsMatrix(:,3) = gui_data.displaySequenceOriginalIndices;
    try
        T = array2table(decisionsMatrix, 'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
        writetable(T, gui_data.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t');
        disp(['Processing decisions saved to: ', gui_data.processingDecisionsFilename]);
    catch ME_write
        rethrow(ME_write);
    end
end