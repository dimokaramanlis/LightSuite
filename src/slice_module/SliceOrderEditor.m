function SliceOrderEditor(optionalVolumePath)
% Interactively reorders slices and adjusts per-slice centering crops.
%
% Keys:
%   Left / Right   navigate slices
%   Enter          move current slice to a new position
%   f              flip current slice left-right
%   o              exclude / restore current slice
%   c              edit the centering crop rectangle (if available)
%   s              save decisions to file
%   Esc            prompt to save and close

    % --- Initial Setup ---
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
        tiffInfo  = imfinfo(inputFileFullPath);
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

    screenSize        = get(0, 'ScreenSize');
    gui_aspect_ratio  = 1.6; gui_width_fraction = 0.5;
    gui_width_px      = screenSize(3) * gui_width_fraction;
    gui_position      = [(screenSize(3)-gui_width_px)/2, ...
                          (screenSize(4)-gui_width_px/gui_aspect_ratio)/2, ...
                          gui_width_px, gui_width_px/gui_aspect_ratio];
    gui_fig = figure('Name', 'Slice Order Editor', 'NumberTitle', 'off', ...
        'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', ...
        'Units', 'pixels', 'Position', gui_position, ...
        'CloseRequestFcn', @(src,evt) callback_close_gui_request(src), ...
        'KeyPressFcn', @callback_keypress);

    gui_data.numSlices                        = numSlices;
    gui_data.currentDisplayPosition           = 1;
    gui_data.flipState                        = zeros(numSlices, 1);
    gui_data.displaySequenceOriginalIndices   = (1:numSlices)';
    gui_data.cropSugg                         = [];   % Nslices x 4: [xmin xmax ymin ymax] in ordering-TIFF pixels
    gui_data.hasCropSugg                      = false;
    gui_data.rectHandle                       = [];
    gui_data.inCropEdit                       = false;

    load_processing_decisions();

    gui_data.imageAxes  = axes('Parent', gui_fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
    gui_data.imageAxes.Colormap = colormap('gray');
    axis(gui_data.imageAxes, 'image', 'off');
    gui_data.imageHandle = image(gui_data.imageAxes, []);
    gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 10);
    gui_data.excludeMarkerHandle = [];
    gui_data.orderTextHandle = text(gui_data.imageAxes, 0, 0, '', ...
        'FontSize', 24, 'Color', 'yellow', 'FontWeight', 'bold', ...
        'BackgroundColor', [0 0 0 0.5], 'VerticalAlignment', 'top');

    guidata(gui_fig, gui_data);
    display_current_slice(gui_fig);
    if gui_data.numSlices > 0
        uiwait(gui_fig);
    elseif ishandle(gui_fig)
        delete(gui_fig);
        disp('No slices to display. GUI closed.');
    end

    % ------------------------------------------------------------------ %
    function load_processing_decisions()
        if ~exist(gui_data.processingDecisionsFilename, 'file') || gui_data.numSlices == 0
            return;
        end
        try
            loadedTable = readtable(gui_data.processingDecisionsFilename, ...
                'Delimiter', '\t', 'ReadVariableNames', true);
            if height(loadedTable) ~= gui_data.numSlices || ...
                    ~all(ismember({'OriginalIndex','FlipState','NewOrderOriginalIndex'}, ...
                    loadedTable.Properties.VariableNames))
                warning('Decision file has mismatching dimensions or columns. Using defaults.');
                return;
            end
            % flip state
            flipData = loadedTable.FlipState;
            if isnumeric(flipData)
                gui_data.flipState = flipData;
            else
                tempFlipState = zeros(gui_data.numSlices, 1);
                for k_f = 1:gui_data.numSlices
                    val_str = strtrim(string(flipData{k_f}));
                    if strcmpi(val_str, "1"),  tempFlipState(k_f) =  1;
                    elseif strcmpi(val_str, "-1"), tempFlipState(k_f) = -1;
                    end
                end
                gui_data.flipState = tempFlipState;
            end
            % ordering
            loadedOrder = loadedTable.NewOrderOriginalIndex;
            if ~any(isnan(loadedOrder)) && numel(unique(loadedOrder)) == gui_data.numSlices
                gui_data.displaySequenceOriginalIndices = loadedOrder;
                disp('Loaded saved order from decisions file.');
            else
                disp('Saved order is incomplete or invalid. Using default 1-N order.');
            end
            % crop suggestions
            cropCols = {'xmin','xmax','ymin','ymax'};
            if all(ismember(cropCols, loadedTable.Properties.VariableNames))
                gui_data.hasCropSugg = true;
                gui_data.cropSugg    = [loadedTable.xmin(:), loadedTable.xmax(:), ...
                                         loadedTable.ymin(:), loadedTable.ymax(:)];
                disp('Loaded centering crop suggestions from decisions file.');
            end
        catch ME_load
            warning('Error loading/parsing decisions file: %s. Using defaults.', ME_load.message);
        end
    end
end

% ======================================================================= %

function display_current_slice(fig)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end

    % clear overlays
    if ~isempty(gui_data.excludeMarkerHandle) && ishandle(gui_data.excludeMarkerHandle)
        delete(gui_data.excludeMarkerHandle); gui_data.excludeMarkerHandle = [];
    end
    if ~isempty(gui_data.rectHandle) && ishandle(gui_data.rectHandle)
        delete(gui_data.rectHandle); gui_data.rectHandle = [];
    end

    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    imgData          = gui_data.originalSliceImages{originalSliceIdx};
    currentFlipState = gui_data.flipState(originalSliceIdx);

    status_str = 'Normal';
    if currentFlipState == 1,  imgData = fliplr(imgData); status_str = 'FLIPPED';
    elseif currentFlipState == -1, status_str = 'EXCLUDED'; end

    set(gui_data.imageHandle, 'CData', imgData);
    axis(gui_data.imageAxes, 'image', 'off');

    if currentFlipState == -1
        [h, w, ~] = size(imgData);
        gui_data.excludeMarkerHandle = text(gui_data.imageAxes, w/2, h/2, 'X', ...
            'FontSize', h/4, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
    end

    % draw crop rectangle
    if gui_data.hasCropSugg && ~isempty(gui_data.cropSugg) && currentFlipState ~= -1
        [imgH, imgW, ~] = size(imgData);
        xmin_c = gui_data.cropSugg(originalSliceIdx, 1);
        xmax_c = gui_data.cropSugg(originalSliceIdx, 2);
        ymin_c = gui_data.cropSugg(originalSliceIdx, 3);
        ymax_c = gui_data.cropSugg(originalSliceIdx, 4);
        if currentFlipState == 1  % image was fliplr — mirror x
            tmp    = imgW - xmax_c + 1;
            xmax_c = imgW - xmin_c + 1;
            xmin_c = tmp;
        end
        xmin_c = max(1, min(imgW, xmin_c));
        xmax_c = max(1, min(imgW, xmax_c));
        ymin_c = max(1, min(imgH, ymin_c));
        ymax_c = max(1, min(imgH, ymax_c));
        if xmax_c > xmin_c && ymax_c > ymin_c
            hold(gui_data.imageAxes, 'on');
            gui_data.rectHandle = rectangle(gui_data.imageAxes, ...
                'Position', [xmin_c, ymin_c, xmax_c-xmin_c, ymax_c-ymin_c], ...
                'EdgeColor', 'y', 'LineWidth', 1.5, 'LineStyle', '--');
            hold(gui_data.imageAxes, 'off');
        end
    end

    orderNumStr = sprintf('%d', gui_data.currentDisplayPosition);
    set(gui_data.orderTextHandle, ...
        'Position', [0.02*size(imgData,2), 0.02*size(imgData,1)], ...
        'String', orderNumStr);

    if gui_data.hasCropSugg
        key_hint = 'Keys: Nav = Left/Right | Reorder = Enter | Crop = c | Save = s';
    else
        key_hint = 'Keys: Nav = Left/Right | Reorder = Enter | Save = s';
    end
    title_str = {sprintf('Slice at Order Position: %d/%d (Original Index: %d) - %s', ...
                         gui_data.currentDisplayPosition, gui_data.numSlices, ...
                         originalSliceIdx, status_str), key_hint};
    set(gui_data.titleHandle, 'String', title_str);
    guidata(fig, gui_data);
end

% ======================================================================= %

function callback_keypress(fig, eventdata)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    if gui_data.inCropEdit; return; end  % block navigation while editing a rectangle

    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    switch eventdata.Key
        case 'leftarrow'
            gui_data.currentDisplayPosition = max(1, gui_data.currentDisplayPosition - 1);
        case 'rightarrow'
            gui_data.currentDisplayPosition = min(gui_data.numSlices, gui_data.currentDisplayPosition + 1);
        case 'f'
            if gui_data.flipState(originalSliceIdx) ~= -1
                gui_data.flipState(originalSliceIdx) = ~gui_data.flipState(originalSliceIdx);
            end
        case 'o'
            if gui_data.flipState(originalSliceIdx) == -1
                gui_data.flipState(originalSliceIdx) = 0;
            else
                gui_data.flipState(originalSliceIdx) = -1;
            end
        case 'c'
            if gui_data.hasCropSugg && ~isempty(gui_data.cropSugg)
                guidata(fig, gui_data);
                edit_crop_rectangle(fig);
                return;
            end
        case {'return', 'enter'}
            reorder_slice_callback(fig); return;
        case 's'
            guidata(fig, gui_data);
            save_processing_decisions(guidata(fig));
            return;
        case 'escape'
            callback_close_gui_request(fig); return;
    end
    guidata(fig, gui_data);
    display_current_slice(fig);
end

% ======================================================================= %

function edit_crop_rectangle(fig)
% Let the user interactively drag/resize the centering crop rectangle.
% The rectangle is shown in display coordinates (fliplr applied when needed).
% Storage is always in the original (unflipped) coordinate system.
    gui_data = guidata(fig);
    orig_idx         = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    currentFlipState = gui_data.flipState(orig_idx);

    if currentFlipState == -1; return; end  % excluded slice — nothing to edit

    imgData = gui_data.originalSliceImages{orig_idx};
    if currentFlipState == 1, imgData = fliplr(imgData); end
    [imgH, imgW, ~] = size(imgData);

    % stored crop in original coords → convert to display coords
    xmin_s = gui_data.cropSugg(orig_idx, 1);
    xmax_s = gui_data.cropSugg(orig_idx, 2);
    ymin_s = gui_data.cropSugg(orig_idx, 3);
    ymax_s = gui_data.cropSugg(orig_idx, 4);
    if currentFlipState == 1
        tmp    = imgW - xmax_s + 1;
        xmax_s = imgW - xmin_s + 1;
        xmin_s = tmp;
    end
    % clip to image bounds before handing to drawrectangle
    xmin_s = max(1, min(imgW-1, xmin_s));
    xmax_s = max(xmin_s+1, min(imgW, xmax_s));
    ymin_s = max(1, min(imgH-1, ymin_s));
    ymax_s = max(ymin_s+1, min(imgH, ymax_s));

    % block other keys, update title
    gui_data.inCropEdit = true;
    if ~isempty(gui_data.rectHandle) && ishandle(gui_data.rectHandle)
        delete(gui_data.rectHandle); gui_data.rectHandle = [];
    end
    set(gui_data.titleHandle, 'String', ...
        {sprintf('Slice %d/%d (Orig. %d) — CROP EDIT', ...
                  gui_data.currentDisplayPosition, gui_data.numSlices, orig_idx), ...
         'Drag / resize the rectangle, then double-click to confirm  (Esc = cancel)'});
    guidata(fig, gui_data);

    % draw interactive rectangle
    roi = drawrectangle(gui_data.imageAxes, ...
        'Position',          [xmin_s, ymin_s, xmax_s-xmin_s, ymax_s-ymin_s], ...
        'Color',             'c', ...
        'LineWidth',         2, ...
        'RotationAngle',     0, ...
        'FixedAspectRatio',  false);
    try
        wait(roi);                       % blocks until user double-clicks
        if isvalid(roi)
            pos = roi.Position;          % [x y width height] in display coords
            delete(roi);

            new_xmin = round(pos(1));
            new_ymin = round(pos(2));
            new_xmax = round(pos(1) + pos(3));
            new_ymax = round(pos(2) + pos(4));

            % clip
            new_xmin = max(1,        min(imgW-1, new_xmin));
            new_xmax = max(new_xmin+1, min(imgW,   new_xmax));
            new_ymin = max(1,        min(imgH-1, new_ymin));
            new_ymax = max(new_ymin+1, min(imgH,   new_ymax));

            % un-flip x if the image was displayed mirrored
            if currentFlipState == 1
                tmp      = imgW - new_xmax + 1;
                new_xmax = imgW - new_xmin + 1;
                new_xmin = tmp;
            end

            gui_data = guidata(fig);
            gui_data.cropSugg(orig_idx, :) = [new_xmin, new_xmax, new_ymin, new_ymax];
            guidata(fig, gui_data);  % persist before the cleanup block re-fetches
        else
            if isvalid(roi), delete(roi); end
        end
    catch
        try; delete(roi); catch; end
    end

    gui_data = guidata(fig);
    gui_data.inCropEdit = false;
    guidata(fig, gui_data);
    display_current_slice(fig);
end

% ======================================================================= %

function reorder_slice_callback(fig)
    gui_data          = guidata(fig);
    current_pos       = gui_data.currentDisplayPosition;
    slice_to_move_idx = gui_data.displaySequenceOriginalIndices(current_pos);

    prompt    = {sprintf('Move slice (Orig. Idx: %d) to new order position (1-%d):', ...
                          slice_to_move_idx, gui_data.numSlices)};
    dlg_title = 'Move to Position';
    answer    = inputdlg(prompt, dlg_title, [1 60], {num2str(current_pos)});
    if isempty(answer), return; end

    new_pos = round(str2double(answer{1}));
    if isnan(new_pos) || new_pos < 1 || new_pos > gui_data.numSlices
        warndlg('Invalid input. Please enter a valid position number.', 'Input Error');
        return;
    end

    new_order_vec = perform_move_to_position( ...
        gui_data.displaySequenceOriginalIndices, slice_to_move_idx, new_pos);
    gui_data.displaySequenceOriginalIndices = new_order_vec;
    gui_data.currentDisplayPosition        = new_pos;

    guidata(fig, gui_data);
    display_current_slice(fig);
end

function new_order_vec = perform_move_to_position(order_vec, slice_to_move_idx, new_pos)
    temp_order_vec = order_vec(order_vec ~= slice_to_move_idx);
    if new_pos == 1
        new_order_vec = [slice_to_move_idx; temp_order_vec];
    elseif new_pos == numel(order_vec)
        new_order_vec = [temp_order_vec; slice_to_move_idx];
    else
        new_order_vec = [temp_order_vec(1:new_pos-1); slice_to_move_idx; temp_order_vec(new_pos:end)];
    end
    disp(sprintf('Moved slice %d to position %d.', slice_to_move_idx, new_pos));
end

% ======================================================================= %

function callback_close_gui_request(fig)
    choice = questdlg('Save changes before closing?', 'Confirm Close', ...
        'Save & Close', 'Discard & Close', 'Cancel', 'Save & Close');
    if strcmp(choice, 'Save & Close'),  save_and_close_confirmed(fig);
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
    decisionsMatrix      = zeros(gui_data.numSlices, 3);
    decisionsMatrix(:,1) = (1:gui_data.numSlices)';
    decisionsMatrix(:,2) = gui_data.flipState;
    decisionsMatrix(:,3) = gui_data.displaySequenceOriginalIndices;
    try
        T = array2table(decisionsMatrix, ...
            'VariableNames', {'OriginalIndex','FlipState','NewOrderOriginalIndex'});
        % persist crop suggestions from memory (may have been edited by user)
        if gui_data.hasCropSugg && ~isempty(gui_data.cropSugg)
            T.xmin = gui_data.cropSugg(:, 1);
            T.xmax = gui_data.cropSugg(:, 2);
            T.ymin = gui_data.cropSugg(:, 3);
            T.ymax = gui_data.cropSugg(:, 4);
        end
        writetable(T, gui_data.processingDecisionsFilename, ...
            'WriteVariableNames', true, 'Delimiter', '\t');
        disp(['Processing decisions saved to: ', gui_data.processingDecisionsFilename]);
    catch ME_write
        rethrow(ME_write);
    end
end
