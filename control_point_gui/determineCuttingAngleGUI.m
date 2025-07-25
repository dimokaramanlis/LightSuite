function determineCuttingAngleGUI(opts)
% determineCuttingAngleGUI allows a user to determine the cutting angle of a
% 3D brain volume by visually matching it to the Allen CCF reference atlas.
%
% This version provides robust view restoration within the GUI while saving
% the correct vector data for post-processing.
%
% --- User-configurable parameters ---
rotation_step = 0.3; % Degrees to rotate per key press
% Initialize the main data structure for the GUI
gui_data = struct;
gui_data.rotation_step = rotation_step;
opts.downfac_reg = opts.allenres/opts.registres;
gui_data.save_path = opts.procpath;
%--------------------------------------------------------------------------
% Load atlas
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
factv       = 255/single(max(gui_data.tv,[],"all"));
gui_data.tv = uint8(single(gui_data.tv(opts.atlasaplims(1):opts.atlasaplims(2), :, :))*factv);
gui_data.tv = imresize3(gui_data.tv,opts.downfac_reg);
gui_data.av = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
gui_data.av = imresize3(gui_data.av(opts.atlasaplims(1):opts.atlasaplims(2), :, :),...
    opts.downfac_reg, "Method","nearest");
gui_data.Rmoving  = imref3d(size(gui_data.av));
disp('Done.')
%--------------------------------------------------------------------------
% Load user's 3D brain volume
volume_dir       = dir(fullfile(opts.procpath,'*inspection.tif*'));
volpath          = fullfile(volume_dir.folder, volume_dir.name);
volload          = readDownStack(volpath);
volload          = permute(volload, [1 2 4 3]);
volload          = permute(volload, [opts.howtoperm 4]);
volload          = single(volload);
minvals          = single(quantile(volload, 0.01, [2 3]));
maxvals          = single(quantile(volload, 0.999, [2 3]));
volload          = 255 * (volload - minvals)./(maxvals - minvals);
gui_data.volume  = uint8(volload);
gui_data.num_user_slices = size(gui_data.volume, 1);
gui_data.colsuse = 1:size(gui_data.volume, 4);
% --- 2. GUI Setup ---
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7;
gui_width_fraction = 0.6;
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [(screen_size_px(3)-gui_width_px)/2, (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, gui_width_px, gui_width_px/gui_aspect_ratio];
gui_fig = figure('KeyPressFcn',@keypress, ...
    'WindowScrollWheelFcn',@scroll_atlas_slice,...
    'Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position, ...
    'CloseRequestFcn',@close_gui);
% Left subplot
gui_data.user_slice_ax = subplot('Position', [0 0 0.5 1]);
hold(gui_data.user_slice_ax, 'on');
axis(gui_data.user_slice_ax, 'image', 'off');
initial_slice_img = squeeze(gui_data.volume(1, :, :, gui_data.colsuse));
gui_data.user_slice_img = image(gui_data.user_slice_ax, initial_slice_img, 'CDataMapping', 'scaled');
colormap(gui_data.user_slice_ax, 'gray');
gui_data.curr_user_slice = 1;
gui_data.user_slice_title = title(gui_data.user_slice_ax, 'Unsaved Slice 1');
gui_data.stats_text = uicontrol('Style', 'text', 'String', 'Marked: 0', ...
    'Units', 'normalized', 'Position', [0.05, 0.88, 0.4, 0.05], ...
    'BackgroundColor', 'w', 'FontSize', 9, 'HorizontalAlignment', 'center');
set(gui_data.user_slice_ax, 'YDir', 'reverse');
% Boundary overlay
gui_data.atlas_boundary_overlay = line(gui_data.user_slice_ax, NaN, NaN, ...
    'Color', 'r', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 4, 'Visible', 'off');
% Right subplot
gui_data.atlas_ax = subplot('Position', [0.55 0 0.45 0.9]);
hold(gui_data.atlas_ax, 'on');
axis(gui_data.atlas_ax, 'vis3d', 'equal', 'manual', 'off');
set(gui_data.atlas_ax, 'Color', 'k', 'ZDir', 'reverse');
view([90, 0]);
[ap_max, dv_max, ml_max] = size(gui_data.tv);
xlim(gui_data.atlas_ax, [1, ap_max]); ylim(gui_data.atlas_ax, [1, ml_max]); zlim(gui_data.atlas_ax, [1, dv_max]);
colormap(gui_data.atlas_ax, 'gray');
clim(gui_data.atlas_ax,[0,255]);
gui_data.atlas_slice_plot = surface(gui_data.atlas_ax, 'EdgeColor', 'none');
% uicontrol for title
gui_data.atlas_title_handle = uicontrol('Style', 'text', 'String', 'Atlas View', ...
    'Units', 'normalized', 'Position', [0.55, 0.92, 0.4, 0.06], ...
    'BackgroundColor', 'w', 'FontSize', 10, 'HorizontalAlignment', 'center');
% Initialize atlas parameters
gui_data.saved_views = repmat({nan(1, 2)}, gui_data.num_user_slices, 1);
gui_data.saved_slice_vectors = repmat({nan(1, 3)}, gui_data.num_user_slices, 1);
gui_data.saved_slice_points = nan(gui_data.num_user_slices, 3);
gui_data.atlas_slice_point = camtarget(gui_data.atlas_ax);

guidata(gui_fig, gui_data);
update_user_slice_view(gui_fig); % Initial view setup
update_stats_display(gui_fig);   % Initial stats display
show_controls_dialog();
end


function keypress(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    is_shift = any(strcmp(eventdata.Modifier, 'shift'));
    step = gui_data.rotation_step;
    switch eventdata.Key
        case 'rightarrow'
            if ~is_shift
                gui_data.curr_user_slice = min(gui_data.curr_user_slice + 1, gui_data.num_user_slices);
                guidata(gui_fig, gui_data);
                update_user_slice_view(gui_fig);
            else
                set(gui_data.atlas_ax, 'View', get(gui_data.atlas_ax, 'View') + [-step, 0]);
                update_atlas_view(gui_fig);
            end
        case 'leftarrow'
            if ~is_shift
                gui_data.curr_user_slice = max(gui_data.curr_user_slice - 1, 1);
                guidata(gui_fig, gui_data);
                update_user_slice_view(gui_fig);
            else
                set(gui_data.atlas_ax, 'View', get(gui_data.atlas_ax, 'View') + [step, 0]);
                update_atlas_view(gui_fig);
            end
        case 'uparrow'
            if is_shift
                set(gui_data.atlas_ax, 'View', get(gui_data.atlas_ax, 'View') + [0, -step]);
                update_atlas_view(gui_fig);
            end
        case 'downarrow'
            if is_shift
                set(gui_data.atlas_ax, 'View', get(gui_data.atlas_ax, 'View') + [0, step]);
                update_atlas_view(gui_fig);
            end
        case 'space'
            current_visibility = get(gui_data.atlas_boundary_overlay, 'Visible');
            if strcmp(current_visibility, 'on')
                set(gui_data.atlas_boundary_overlay, 'Visible', 'off');
            else
                set(gui_data.atlas_boundary_overlay, 'Visible', 'on');
                update_boundary_overlay(gui_fig);
            end
        case 'return'
            % Save the current view and position for this slice
            gui_data.saved_views{gui_data.curr_user_slice} = get(gui_data.atlas_ax, 'View');
            gui_data.saved_slice_vectors{gui_data.curr_user_slice} = get_camera_vector(gui_data);
            gui_data.saved_slice_points(gui_data.curr_user_slice, :) = gui_data.atlas_slice_point;
            guidata(gui_fig, gui_data);
            update_user_slice_view(gui_fig);
            update_stats_display(gui_fig);

        %% MODIFICATION START: Add a case to clear the current slice's saved data
        case 'c'
            slice_idx = gui_data.curr_user_slice;
            % Reset data for the current slice to default 'unsaved' state
            gui_data.saved_views{slice_idx} = nan(1, 2);
            gui_data.saved_slice_vectors{slice_idx} = nan(1, 3);
            gui_data.saved_slice_points(slice_idx, :) = nan(1, 3);
            
            disp(['Cleared saved position for slice ' num2str(slice_idx)]);
            guidata(gui_fig, gui_data);
            
            % Refresh the view and stats display
            update_user_slice_view(gui_fig);
            update_stats_display(gui_fig);
        %% MODIFICATION END

        case {'1', '2', '3'}
             if str2double(eventdata.Key) <= size(gui_data.volume, 4)
                gui_data.colsuse = str2double(eventdata.Key);
                update_atlas_view(gui_fig);
                guidata(gui_fig,gui_data);
                update_user_slice_view(gui_fig);
             end
        case '0'
            gui_data.colsuse = 1:size(gui_data.volume, 4);
            update_atlas_view(gui_fig);
            guidata(gui_fig,gui_data);
            update_user_slice_view(gui_fig);
    end
end
function scroll_atlas_slice(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    cam_vector = get_camera_vector(gui_data);
    scroll_speed_factor = 3;
    gui_data.atlas_slice_point = gui_data.atlas_slice_point + (eventdata.VerticalScrollCount * cam_vector * scroll_speed_factor);
    guidata(gui_fig, gui_data);
    update_atlas_view(gui_fig);
end
function close_gui(gui_fig, ~)
    gui_data = guidata(gui_fig);
    if any(isnan(gui_data.saved_slice_points(:)))
        uiwait(msgbox('Warning: Not all slices have a saved atlas position.', 'Incomplete', 'warn'));
    end
    answer = questdlg('Save the alignment data?', 'Confirm Exit', 'Yes', 'No', 'Cancel', 'Yes');
    if strcmp(answer, 'Yes')
        disp('Saving alignment data...');
        cutting_angle_data = struct();
        cutting_angle_data.saved_slice_vectors = gui_data.saved_slice_vectors;
        cutting_angle_data.slice_points = gui_data.saved_slice_points;
        %% MODIFICATION START: Save view data to allow for full session restoration
        cutting_angle_data.saved_views = gui_data.saved_views;
        %% MODIFICATION END
        save_filename = fullfile(gui_data.save_path, 'cutting_angle_data.mat');
        save(save_filename, 'cutting_angle_data');
        fprintf('Data saved to: %s\n', save_filename);
        delete(gui_fig);
    elseif strcmp(answer, 'No')
        disp('Closing without saving.');
        delete(gui_fig);
    end
end

%% MODIFICATION START: Add helper function to get the mean of saved view angles
function mean_view = get_mean_view(gui_data)
    % Find which views have been saved
    is_saved = ~cellfun(@(v) all(isnan(v)), gui_data.saved_views);
    
    if ~any(is_saved)
        % If no slices are saved yet, return a default view
        %ANAS - MODIFICATION: if nothing is saved, keep the current view.
        mean_view = get(gui_data.atlas_ax, 'View');%[90, 0]; 
        return;
    end
    
    % Get the matrix of saved [azimuth, elevation] pairs
    saved_views_matrix = cell2mat(gui_data.saved_views(is_saved));
    
    % Calculate the mean view
    mean_view = mean(saved_views_matrix, 1);
end
%% MODIFICATION END

function update_stats_display(gui_fig)
    gui_data = guidata(gui_fig);
    is_saved = ~cellfun(@(v) all(isnan(v)), gui_data.saved_views);
    num_marked = sum(is_saved);
    if num_marked == 0
        stats_str = sprintf('Marked: 0 / %d', gui_data.num_user_slices);
    else
        saved_views_matrix = cell2mat(gui_data.saved_views(is_saved));
        dv_angles = saved_views_matrix(:, 1) - 90;
        ml_angles = saved_views_matrix(:, 2);
        avg_dv_angle = mean(dv_angles);
        avg_ml_angle = mean(ml_angles);
        stats_str = sprintf('Marked: %d / %d | Avg Angle (DV: %.1f°, ML: %.1f°)', ...
                          num_marked, gui_data.num_user_slices, avg_dv_angle, avg_ml_angle);
    end
    set(gui_data.stats_text, 'String', stats_str);
end

%% MODIFICATION START: Update function to use learning behavior for unsaved slices
function update_user_slice_view(gui_fig)
    gui_data = guidata(gui_fig);
    slice_idx = gui_data.curr_user_slice;
    set(gui_data.user_slice_img, 'CData', squeeze(gui_data.volume(slice_idx, :, :, gui_data.colsuse)));

    % Check if the current slice has a saved position
    if all(~isnan(gui_data.saved_slice_points(slice_idx, :)))
        % --- SLICE IS SAVED ---
        % Restore the specific saved view angles
        saved_view = gui_data.saved_views{slice_idx};
        view(gui_data.atlas_ax, saved_view);
        
        % Restore the saved slice plane position
        gui_data.atlas_slice_point = gui_data.saved_slice_points(slice_idx, :);
        guidata(gui_fig, gui_data);
        
        % Update title to show "SAVED"
        gui_data.user_slice_title.String = sprintf('Slice %d - SAVED', slice_idx);
        gui_data.user_slice_title.Color = 'g';
        
    else
        % --- SLICE IS UNSAVED (APPLY LEARNING) ---
        % Calculate the mean angle from all other saved slices
        mean_view = get_mean_view(gui_data);
        
        % Set the atlas view to the predicted mean angle
        view(gui_data.atlas_ax, mean_view);
        
        % Update title to show "Unsaved"
        gui_data.user_slice_title.String = sprintf('Unsaved Slice %d', slice_idx);
        gui_data.user_slice_title.Color = 'k';
    end

    % Update the atlas slice plane and boundary overlay for both cases
    update_atlas_view(gui_fig);
    update_boundary_overlay(gui_fig);
end
%% MODIFICATION END

function update_atlas_view(gui_fig)
    gui_data = guidata(gui_fig);
    
    [tv_slice, ~, plane_coords] = get_atlas_slice(gui_data, 3);
    
    set(gui_data.atlas_slice_plot, 'XData', plane_coords.ap, 'YData', plane_coords.ml, 'ZData', plane_coords.dv, 'CData', tv_slice);
        
    % Update atlas title
    [az, el] = view(gui_data.atlas_ax);
    angle_dv_axis = az - 90;
    angle_ml_axis = el;
    title_str1 = sprintf('Angle (DV-axis): %.1f°, Angle (ML-axis): %.1f°', angle_dv_axis, angle_ml_axis);
    
    cam_vector = get_camera_vector(gui_data);
    plane_offset = -(cam_vector * gui_data.atlas_slice_point');
    plane_offset_mm = plane_offset / 100;
    title_str2 = sprintf('Slice Position: %.2f mm', plane_offset_mm);
    
    set(gui_data.atlas_title_handle, 'String', {title_str1, title_str2});
end
function update_boundary_overlay(gui_fig)
    gui_data = guidata(gui_fig);
    if strcmp(get(gui_data.atlas_boundary_overlay, 'Visible'), 'off')
        set(gui_data.atlas_boundary_overlay, 'XData', NaN, 'YData', NaN);
        return;
    end
    
    [~, av_slice] = get_atlas_slice(gui_data, 1);
    
    if all(isnan(av_slice), 'all')
        set(gui_data.atlas_boundary_overlay, 'XData', NaN, 'YData', NaN);
        return;
    end

    av_boundaries = round(conv2(av_slice, ones(2)./4, 'same')) ~= av_slice;
    [rows, cols] = find(av_boundaries);
    
    [user_h, user_w] = size(squeeze(gui_data.volume(gui_data.curr_user_slice, :, :, 1)));
    [atlas_h, atlas_w] = size(av_slice);
    scale_h = user_h / atlas_h;
    scale_w = user_w / atlas_w;
    
    set(gui_data.atlas_boundary_overlay, 'XData', cols * scale_w, 'YData', rows * scale_h);
end
function cam_vector = get_camera_vector(gui_data)
    campos_val = campos(gui_data.atlas_ax);
    camtarget_val = camtarget(gui_data.atlas_ax);
    cam_vector = (camtarget_val - campos_val) ./ norm(camtarget_val - campos_val);
end
function [tv_slice, av_slice, plane_coords] = get_atlas_slice(gui_data, spacing)
    cam_vector = get_camera_vector(gui_data);
    plane_offset = -(cam_vector * gui_data.atlas_slice_point');
    [~, main_axis] = max(abs(cam_vector));
    s = size(gui_data.tv);
    
    switch main_axis
        case 1
            [ml_grid, dv_grid] = meshgrid(1:spacing:s(3), 1:spacing:s(2));
            ap_grid = (-plane_offset - cam_vector(2)*ml_grid - cam_vector(3)*dv_grid) / cam_vector(1);
        case 2
            [ap_grid, dv_grid] = meshgrid(1:spacing:s(1), 1:spacing:s(2));
            ml_grid = (-plane_offset - cam_vector(1)*ap_grid - cam_vector(3)*dv_grid) / cam_vector(2);
        case 3
            [ap_grid, ml_grid] = meshgrid(1:spacing:s(1), 1:spacing:s(3));
            dv_grid = (-plane_offset - cam_vector(1)*ap_grid - cam_vector(2)*ml_grid) / cam_vector(3);
    end
    plane_coords = struct('ap', ap_grid, 'ml', ml_grid, 'dv', dv_grid);
    
    ap_idx = round(ap_grid); ml_idx = round(ml_grid); dv_idx = round(dv_grid);
    valid_mask = ap_idx > 0 & ap_idx <= s(1) & dv_idx > 0 & dv_idx <= s(2) & ml_idx > 0 & ml_idx <= s(3);
    indices = sub2ind(s, ap_idx(valid_mask), dv_idx(valid_mask), ml_idx(valid_mask));
    
    tv_slice = nan(size(ap_grid));
    av_slice = nan(size(ap_grid));
    tv_slice(valid_mask) = gui_data.tv(indices);
    av_slice(valid_mask) = gui_data.av(indices);
end

%% MODIFICATION START: Update controls dialog to include the new 'c' key
function show_controls_dialog()
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'non-modal';
    msgbox( ...
        {'\fontsize{12}' ...
        '\bf Controls: \rm' ...
        'Left/Right Arrows: Change brain slice', ...
        'Shift + Arrows: Rotate 3D atlas', ...
        'Spacebar: Toggle atlas boundary overlay', ...
        'Scroll Wheel: Move atlas slice plane', ...
        '1-3 or 0: Toggle channels (0 for all)',...
        'Enter: Save atlas position for current slice', ...
        '''c'': Clear saved position for current slice'}, ...
        'Controls',CreateStruct);
end
%% MODIFICATION END