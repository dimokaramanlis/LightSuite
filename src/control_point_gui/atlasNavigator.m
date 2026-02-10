function atlasNavigator
% A simple GUI to navigate through the Allen Brain Atlas using mouse scroll.
% Press the spacebar to toggle annotation outlines.

% --- GUI data structure ---
gui_data = struct;

% --- Load Atlas Data ---
% Ensure the Allen CCF atlas files are in your MATLAB path.
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
if isempty(allen_atlas_path)
    error('Allen CCF atlas not found. Please add the atlas to your MATLAB path.');
end

disp('Loading Allen CCF atlas...');
% Load the average template and annotation volumes
% Note: permute is used to make the AP axis the first dimension for consistency
% with the original code's navigation logic.
gui_data.tv = niftiread(fullfile(allen_atlas_path, 'average_template_10.nii.gz'));
gui_data.tv = imresize3(gui_data.tv, 0.5);
gui_data.av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
gui_data.av = imresize3(gui_data.av, 0.5, 'nearest');
gui_data.tv = uint8(255*single(gui_data.tv)/single(quantile(gui_data.tv, 0.99, 'all')));
disp('Atlas loaded.');

% --- Initialize GUI State ---
gui_data.current_slice = 1;
gui_data.max_slice = size(gui_data.tv, 1);
gui_data.show_annotations = false; % Annotations are initially off

% --- Create GUI Components ---
% Create the main figure and set its callback functions
gui_data.fig = figure('Name', 'Atlas Navigator', 'NumberTitle', 'off', ...
    'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', ...
    'WindowScrollWheelFcn', @scroll_slice_callback, ...
    'KeyPressFcn', @key_press_callback);

% Create an axes for displaying the atlas slices
gui_data.atlas_ax = axes('Parent', gui_data.fig, 'Position', [0.05, 0.05, 0.9, 0.9]);
gui_data.atlas_ax.YDir = 'reverse';
axis image off;
hold(gui_data.atlas_ax, 'on');
colormap(gui_data.atlas_ax, 'gray');

% Display the first slice of the atlas
initial_slice_img = squeeze(gui_data.tv(gui_data.current_slice, :, :));
gui_data.atlas_im = image(gui_data.atlas_ax, initial_slice_img);

% Create an empty plot handle for the annotation outlines
gui_data.annotation_outline = plot(gui_data.atlas_ax, NaN, NaN, 'r.', 'MarkerSize', 2);

% Add a title that will be updated
gui_data.atlas_title = title(gui_data.atlas_ax, '', 'FontSize', 14);

% Store the GUI data and perform the initial update
guidata(gui_data.fig, gui_data);
update_display(gui_data.fig);

%% --- Callback Functions ---

    function scroll_slice_callback(hObject, eventdata)
        % This function is called when the mouse wheel is scrolled
        
        % Get the current GUI data
        gui_data = guidata(hObject);
        
        % Get the scroll direction and update the slice number
        scroll_count = eventdata.VerticalScrollCount;
        gui_data.current_slice = gui_data.current_slice + 5 * scroll_count;
        
        % Clamp the slice number to be within the valid range
        gui_data.current_slice = max(1, gui_data.current_slice);
        gui_data.current_slice = min(gui_data.max_slice, gui_data.current_slice);
        
        % Store the updated GUI data and refresh the display
        guidata(hObject, gui_data);
        update_display(hObject);
    end

    function key_press_callback(hObject, eventdata)
        % This function is called when a key is pressed
        
        % Get the current GUI data
        gui_data = guidata(hObject);
        
        % Toggle annotation visibility if the spacebar is pressed
        if strcmp(eventdata.Key, 'space')
            gui_data.show_annotations = ~gui_data.show_annotations;
        end
        
        % Store the updated GUI data and refresh the display
        guidata(hObject, gui_data);
        update_display(hObject);
    end

%% --- Helper Functions ---

    function update_display(fig_handle)
        % This function refreshes the entire GUI display
        
        % Get the current GUI data
        gui_data = guidata(fig_handle);
        
        % Round the current slice to ensure it's an integer index
        slice_idx = round(gui_data.current_slice);
        
        % Update the main atlas image
        set(gui_data.atlas_im, 'CData', squeeze(gui_data.tv(slice_idx, :, :)));
        
        % Update the title
        set(gui_data.atlas_title, 'String', ...
            sprintf('AP Slice: %d / %d. Press SPACE for annotations', slice_idx, gui_data.max_slice));
        
        % Update the annotation outlines based on visibility state
        if gui_data.show_annotations
            % Get the annotation slice
            av_slice = squeeze(gui_data.av(slice_idx, :, :));
            
            % Find the boundaries of all annotated regions
            % We use edge detection on the label matrix to find where regions meet
            boundaries = edge(av_slice, 'canny');
            
            % Get the coordinates of the boundary pixels
            [rows, cols] = find(boundaries);
            
            % Update the plot with the new boundary coordinates
            % We plot (cols, rows) because imagesc displays matrix C(m,n) at (x,y) = (n,m)
            set(gui_data.annotation_outline, 'XData', cols, 'YData', rows);
        else
            % If annotations are hidden, clear the plot data
            set(gui_data.annotation_outline, 'XData', NaN, 'YData', NaN);
        end
    end

end