function SliceOrderEditor(optionalVolumePath)
% Interactively reorders slices and adjusts per-slice centering crops.
%
% Keys:
%   Left / Right   navigate slices
%   Enter          move current slice to a new position
%   f              flip current slice left-right
%   o              exclude / restore current slice
%   c              edit the centering crop rectangle (if available)
%   m              show / hide the montage of all slices
%   s              save decisions to file
%   Esc            prompt to save and close
%
% In the montage, click a tile to show that slice, or drag it onto another
% tile to move the slice there.

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
    % montage window (its code is at the end of this file)
    gui_data.showMontage                      = true;   % open the montage at start
    gui_data.montageColumns                   = 8;
    gui_data.montageThumbWidth                = 150;    % tile width in pixels
    gui_data.montageFig                       = [];
    gui_data.montageAxes                      = [];
    gui_data.montageImage                     = [];
    gui_data.montageHighlight                 = [];
    gui_data.montageTarget                    = [];     % drop target while dragging
    gui_data.montageLabels                    = [];
    gui_data.montageCrops                     = [];     % crop rectangle on each tile
    gui_data.montageMarks                     = [];     % red crosses on excluded slices
    gui_data.montageThumbs                    = {};
    gui_data.dragFromPosition                 = [];

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
    build_montage_window(gui_fig);   % montage
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
        key_hint = 'Keys: Nav = Left/Right | Reorder = Enter | Crop = c | Montage = m | Save = s';
    else
        key_hint = 'Keys: Nav = Left/Right | Reorder = Enter | Montage = m | Save = s';
    end
    title_str = {sprintf('Slice at Order Position: %d/%d (Original Index: %d) - %s', ...
                         gui_data.currentDisplayPosition, gui_data.numSlices, ...
                         originalSliceIdx, status_str), key_hint};
    set(gui_data.titleHandle, 'String', title_str);
    guidata(fig, gui_data);
    update_montage_highlight(fig);   % montage
end

% ======================================================================= %

function callback_keypress(fig, eventdata)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    if gui_data.inCropEdit; return; end  % block navigation while editing a rectangle

    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    refresh_tiles = false;   % montage: f and o change the tile
    switch eventdata.Key
        case 'leftarrow'
            gui_data.currentDisplayPosition = max(1, gui_data.currentDisplayPosition - 1);
        case 'rightarrow'
            gui_data.currentDisplayPosition = min(gui_data.numSlices, gui_data.currentDisplayPosition + 1);
        case 'f'
            if gui_data.flipState(originalSliceIdx) ~= -1
                gui_data.flipState(originalSliceIdx) = ~gui_data.flipState(originalSliceIdx);
            end
            refresh_tiles = true;
        case 'o'
            if gui_data.flipState(originalSliceIdx) == -1
                gui_data.flipState(originalSliceIdx) = 0;
            else
                gui_data.flipState(originalSliceIdx) = -1;
            end
            refresh_tiles = true;
        case 'm'
            toggle_montage_window(fig); return;   % montage
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
    if refresh_tiles, refresh_montage(fig); end   % montage
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
    refresh_montage(fig);   % montage
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
    refresh_montage(fig);   % montage
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
    elseif strcmp(choice, 'Discard & Close')
        close_montage_window(fig);   % montage
        if ishandle(fig); delete(fig); end
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
    close_montage_window(fig);   % montage
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
    save_montage_snapshot(gui_data);   % montage
end

% ======================================================================= %
%  Montage window: all slices in the current order, one tile each
% ======================================================================= %

function build_montage_window(fig)
    gui_data = guidata(fig);
    if ~gui_data.showMontage || gui_data.numSlices == 0, return; end

    % thumbnails are made once and kept when the window is reopened
    if isempty(gui_data.montageThumbs)
        gui_data.montageThumbs = make_thumbnails(gui_data.originalSliceImages, gui_data.montageThumbWidth);
    end

    [th, tw, ~] = size(gui_data.montageThumbs{1});
    n_cols = gui_data.montageColumns;
    n_rows = ceil(gui_data.numSlices / n_cols);

    % size the window to the montage, but keep it on screen
    screen = get(0, 'ScreenSize');
    win_w  = min(n_cols * tw + 40, screen(3) * 0.9);
    win_h  = min(n_rows * th + 60, screen(4) * 0.85);

    gui_data.montageFig = figure( ...
        'Name', 'Slice Order Editor - Montage', 'NumberTitle', 'off', ...
        'Toolbar', 'none', 'Menubar', 'none', 'Color', 'k', 'Units', 'pixels', ...
        'Position', [40, max(60, screen(4) - win_h - 100), win_w, win_h], ...
        'CloseRequestFcn', @(src, evt) montage_close_request(src), ...
        'KeyPressFcn', @(src, evt) forward_key_to_main(fig, evt), ...
        'WindowButtonDownFcn', @(src, evt) montage_button_down(fig));

    gui_data.montageAxes = axes('Parent', gui_data.montageFig, ...
        'Units', 'normalized', 'Position', [0 0 1 1], 'Color', 'k');
    gui_data.montageImage = image(gui_data.montageAxes, zeros(n_rows*th, n_cols*tw, 3, 'uint8'));
    axis(gui_data.montageAxes, 'image', 'off');
    hold(gui_data.montageAxes, 'on');

    % frame on the slice shown in the main window
    gui_data.montageHighlight = rectangle('Parent', gui_data.montageAxes, ...
        'Position', [0.5 0.5 tw th], 'EdgeColor', [1 1 0], 'LineWidth', 3);

    % drop target, shown only while dragging
    gui_data.montageTarget = rectangle('Parent', gui_data.montageAxes, ...
        'Position', [0.5 0.5 tw th], 'EdgeColor', [0 1 1], 'LineWidth', 2, ...
        'LineStyle', '--', 'Visible', 'off');

    % one label and one crop rectangle per tile, filled in by refresh_montage
    gui_data.montageLabels = gobjects(gui_data.numSlices, 1);
    gui_data.montageCrops  = gobjects(gui_data.numSlices, 1);
    for k = 1:gui_data.numSlices
        [x0, y0] = tile_corner(k, n_cols, th, tw);
        gui_data.montageCrops(k) = rectangle('Parent', gui_data.montageAxes, ...
            'Position', [x0 y0 1 1], 'EdgeColor', 'y', 'LineWidth', 1, ...
            'LineStyle', '--', 'Visible', 'off', 'HitTest', 'off');
        gui_data.montageLabels(k) = text(gui_data.montageAxes, x0 + 4, y0 + 10, '', ...
            'Color', [1 1 0], 'FontSize', 9, 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.5], ...   % readable over the crop line
            'VerticalAlignment', 'top', 'Interpreter', 'none', 'HitTest', 'off');
    end

    guidata(fig, gui_data);
    refresh_montage(fig);

    % give the keyboard back to the main window
    if ishandle(fig), figure(fig); end
end

% ======================================================================= %

function refresh_montage(fig)
% Redraw all tiles from the current order, flips, exclusions and crops.
    gui_data = guidata(fig);
    if ~gui_data.showMontage || isempty(gui_data.montageFig) || ~ishandle(gui_data.montageFig)
        return;
    end

    [th, tw, ~] = size(gui_data.montageThumbs{1});
    n_cols = gui_data.montageColumns;
    n_rows = ceil(gui_data.numSlices / n_cols);

    canvas = zeros(n_rows * th, n_cols * tw, 3, 'uint8');

    % clear the previous red crosses
    delete(gui_data.montageMarks(ishandle(gui_data.montageMarks)));
    gui_data.montageMarks = gobjects(0);

    for pos = 1:gui_data.numSlices
        orig_idx = gui_data.displaySequenceOriginalIndices(pos);
        thumb    = gui_data.montageThumbs{orig_idx};
        state    = gui_data.flipState(orig_idx);

        if state == 1
            thumb = fliplr(thumb);
        elseif state == -1
            thumb = uint8(single(thumb) * 0.35);   % dim excluded slices
        end

        [x0, y0] = tile_corner(pos, n_cols, th, tw);
        canvas(y0:(y0+th-1), x0:(x0+tw-1), :) = thumb;

        if state == -1
            gui_data.montageMarks(end+1) = text(gui_data.montageAxes, ...
                x0 + tw/2, y0 + th/2, 'X', 'Color', [1 0.2 0.2], 'FontSize', 22, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'HitTest', 'off'); %#ok<AGROW>
        end

        label = sprintf('%d [%d]', pos, orig_idx);
        if state == 1, label = [label ' F']; end
        set(gui_data.montageLabels(pos), 'String', label);

        % the slice's crop, as in the main window (none on excluded slices)
        rect = [];
        if gui_data.hasCropSugg && ~isempty(gui_data.cropSugg) && state ~= -1
            rect = tile_crop_rectangle(gui_data, orig_idx, x0, y0);
        end
        if isempty(rect)
            set(gui_data.montageCrops(pos), 'Visible', 'off');
        else
            set(gui_data.montageCrops(pos), 'Position', rect, 'Visible', 'on');
        end
    end

    set(gui_data.montageImage, 'CData', canvas);
    guidata(fig, gui_data);
    update_montage_highlight(fig);
end

% ======================================================================= %

function update_montage_highlight(fig)
% Move the yellow frame to the slice shown in the main window.
    gui_data = guidata(fig);
    if ~gui_data.showMontage || isempty(gui_data.montageHighlight) || ~ishandle(gui_data.montageHighlight)
        return;
    end

    [th, tw, ~] = size(gui_data.montageThumbs{1});
    [x0, y0] = tile_corner(gui_data.currentDisplayPosition, gui_data.montageColumns, th, tw);
    set(gui_data.montageHighlight, 'Position', [x0 - 0.5, y0 - 0.5, tw, th]);
end

% ======================================================================= %

function montage_button_down(fig)
% A press starts a drag; montage_release decides if it was a click or a move.
    gui_data = guidata(fig);
    if gui_data.numSlices == 0, return; end
    if gui_data.inCropEdit, return; end   % not while a crop is being edited

    pos = position_under_cursor(gui_data);
    if isempty(pos), return; end

    gui_data.dragFromPosition = pos;
    guidata(fig, gui_data);

    set(gui_data.montageFig, ...
        'WindowButtonMotionFcn', @(src, evt) montage_drag(fig), ...
        'WindowButtonUpFcn',     @(src, evt) montage_release(fig));
end

function montage_drag(fig)
    gui_data = guidata(fig);
    if isempty(gui_data.dragFromPosition), return; end

    pos = position_under_cursor(gui_data);
    if isempty(pos) || pos == gui_data.dragFromPosition
        set(gui_data.montageTarget, 'Visible', 'off');
        return;
    end

    [th, tw, ~] = size(gui_data.montageThumbs{1});
    [x0, y0] = tile_corner(pos, gui_data.montageColumns, th, tw);
    set(gui_data.montageTarget, 'Position', [x0 - 0.5, y0 - 0.5, tw, th], 'Visible', 'on');
end

function montage_release(fig)
    gui_data = guidata(fig);
    from_pos = gui_data.dragFromPosition;
    if isempty(from_pos), return; end

    set(gui_data.montageFig, 'WindowButtonMotionFcn', '', 'WindowButtonUpFcn', '');
    set(gui_data.montageTarget, 'Visible', 'off');
    gui_data.dragFromPosition = [];
    guidata(fig, gui_data);

    to_pos = position_under_cursor(gui_data);

    if isempty(to_pos) || to_pos == from_pos
        % released on the same tile: a click, show that slice
        gui_data.currentDisplayPosition = from_pos;
        guidata(fig, gui_data);
        display_current_slice(fig);
        return;
    end

    % released on another tile: move the slice there
    slice_to_move = gui_data.displaySequenceOriginalIndices(from_pos);
    gui_data.displaySequenceOriginalIndices = ...
        perform_move_to_position(gui_data.displaySequenceOriginalIndices, slice_to_move, to_pos);
    gui_data.currentDisplayPosition = to_pos;
    guidata(fig, gui_data);

    display_current_slice(fig);
    refresh_montage(fig);
end

% ======================================================================= %

function pos = position_under_cursor(gui_data)
% Order position under the mouse, empty if it is off the grid.
    pos = [];
    if isempty(gui_data.montageAxes) || ~ishandle(gui_data.montageAxes), return; end

    cp = get(gui_data.montageAxes, 'CurrentPoint');
    x  = cp(1, 1);
    y  = cp(1, 2);

    [th, tw, ~] = size(gui_data.montageThumbs{1});
    col = floor((x - 0.5) / tw);
    row = floor((y - 0.5) / th);
    if col < 0 || col >= gui_data.montageColumns || row < 0, return; end

    candidate = row * gui_data.montageColumns + col + 1;
    if candidate >= 1 && candidate <= gui_data.numSlices
        pos = candidate;
    end
end

function [x0, y0] = tile_corner(pos, n_cols, th, tw)
% Top-left pixel of the tile at order position pos.
    row = floor((pos - 1) / n_cols);
    col = mod(pos - 1, n_cols);
    x0  = col * tw + 1;
    y0  = row * th + 1;
end

function rect = tile_crop_rectangle(gui_data, orig_idx, x0, y0)
% The slice's crop scaled to its tile, mirrored if the slice is flipped (as in
% display_current_slice). Empty if the crop has no area.
    [imgH, imgW, ~] = size(gui_data.originalSliceImages{orig_idx});
    xmin_c = gui_data.cropSugg(orig_idx, 1);
    xmax_c = gui_data.cropSugg(orig_idx, 2);
    ymin_c = gui_data.cropSugg(orig_idx, 3);
    ymax_c = gui_data.cropSugg(orig_idx, 4);
    if gui_data.flipState(orig_idx) == 1
        tmp    = imgW - xmax_c + 1;
        xmax_c = imgW - xmin_c + 1;
        xmin_c = tmp;
    end
    xmin_c = max(1, min(imgW, xmin_c));
    xmax_c = max(1, min(imgW, xmax_c));
    ymin_c = max(1, min(imgH, ymin_c));
    ymax_c = max(1, min(imgH, ymax_c));

    rect = [];
    if xmax_c <= xmin_c || ymax_c <= ymin_c, return; end

    scale = gui_data.montageThumbWidth / imgW;   % as in make_thumbnails
    rect  = [x0 - 0.5 + (xmin_c - 0.5) * scale, y0 - 0.5 + (ymin_c - 0.5) * scale, ...
             (xmax_c - xmin_c) * scale, (ymax_c - ymin_c) * scale];
end

function thumbs = make_thumbnails(images, target_width)
% Resize every slice to the same width and pad them to the same height.
    n      = numel(images);
    scaled = cell(n, 1);
    for k = 1:n
        img = images{k};
        if size(img, 3) == 1
            img = repmat(img, 1, 1, 3);   % grayscale stacks
        end
        scale = target_width / size(img, 2);
        scaled{k} = imresize(img, scale);
    end

    max_h  = max(cellfun(@(x) size(x, 1), scaled));
    thumbs = cell(n, 1);
    for k = 1:n
        pad       = max_h - size(scaled{k}, 1);
        thumbs{k} = padarray(scaled{k}, [pad 0 0], 0, 'post');
    end
end

% ======================================================================= %

function forward_key_to_main(fig, evt)
% Keys pressed in the montage go to the main window's handler.
    if ishandle(fig)
        callback_keypress(fig, evt);
    end
end

function montage_close_request(montage_fig)
% Closing the montage leaves the editor open; m reopens it.
    if ishandle(montage_fig), delete(montage_fig); end
end

function save_montage_snapshot(gui_data)
% Save a picture of the montage next to the decisions file.
    if ~isfield(gui_data, 'montageFig') || isempty(gui_data.montageFig) || ~ishandle(gui_data.montageFig)
        disp('Montage window is closed, snapshot not saved (press m to reopen it).');
        return;
    end

    out_dir       = fileparts(gui_data.processingDecisionsFilename);
    snapshot_file = fullfile(out_dir, 'slice_order_montage.png');
    try
        exportgraphics(gui_data.montageFig, snapshot_file, ...
            'Resolution', 120, 'BackgroundColor', 'k');
        disp(['Montage snapshot saved to: ', snapshot_file]);
    catch ME_snap
        % the decisions are already saved, so only warn
        warning('Could not save montage snapshot: %s', ME_snap.message);
    end
end

function toggle_montage_window(fig)
% m: close the montage if it is open, open it otherwise.
    gui_data = guidata(fig);
    if ~isempty(gui_data.montageFig) && ishandle(gui_data.montageFig)
        delete(gui_data.montageFig);
        gui_data.montageFig = [];
        guidata(fig, gui_data);
    else
        gui_data.showMontage = true;   % also when it was off at start
        guidata(fig, gui_data);
        build_montage_window(fig);
    end
end

function close_montage_window(fig)
    if ~ishandle(fig), return; end
    gui_data = guidata(fig);
    if isfield(gui_data, 'montageFig') && ~isempty(gui_data.montageFig) && ishandle(gui_data.montageFig)
        delete(gui_data.montageFig);
    end
end
