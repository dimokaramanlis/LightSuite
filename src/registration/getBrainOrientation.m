function permvec = getBrainOrientation(backvol, atlasvol)
%GETBRAINORIENTATION Interactive GUI to determine the volume permutation
%   permvec = getBrainOrientation(backvol, atlasvol) opens an interactive figure 
%   allowing the user to map the sample's axes to the target atlas.
%   The top row shows the fixed atlas reference.
%   The bottom row dynamically updates the sample based on the selected mapping.

    % Default output if user closes without saving
    permvec = [1 2 3];
    
    % Configuration for dropdown options
    opts_str = {'Dim 1 (+)', 'Dim 1 (Flipped -)', ...
                'Dim 2 (+)', 'Dim 2 (Flipped -)', ...
                'Dim 3 (+)', 'Dim 3 (Flipped -)'};
    opts_val = [1, -1, 2, -2, 3, -3];

    % Create the main figure
    fig = figure('Name', 'Match Sample Orientation to Atlas', ...
                 'Position', [100 100 1400 750], ...
                 'MenuBar', 'none', 'ToolBar', 'none', ...
                 'Color', 'w', 'NumberTitle', 'off', ...
                 'CloseRequestFcn', @closeFigure);

    % --- 1. Bottom Control Panel ---
    controlPanel = uipanel(fig, 'Position', [0.0 0.0 1.0 0.20], ...
                           'Title', 'Orientation Mapping (Match Sample to Atlas)', ...
                           'BackgroundColor', 'w', 'FontName', 'Arial', 'FontWeight', 'bold');

    % Dropdown for Atlas Dim 1
    uicontrol(controlPanel, 'Style', 'text', 'String', 'Map to Atlas Dim 1 (Y-axis):', ...
              'Units', 'normalized', 'Position', [0.02 0.6 0.2 0.2], ...
              'HorizontalAlignment', 'left', 'FontSize', 11, 'BackgroundColor', 'w');
    dd1 = uicontrol(controlPanel, 'Style', 'popupmenu', 'String', opts_str, ...
                    'Units', 'normalized', 'Position', [0.22 0.6 0.15 0.2], 'Value', 1, ...
                    'Callback', @updatePreview, 'FontSize', 11);

    % Dropdown for Atlas Dim 2
    uicontrol(controlPanel, 'Style', 'text', 'String', 'Map to Atlas Dim 2 (X-axis):', ...
              'Units', 'normalized', 'Position', [0.02 0.35 0.2 0.2], ...
              'HorizontalAlignment', 'left', 'FontSize', 11, 'BackgroundColor', 'w');
    dd2 = uicontrol(controlPanel, 'Style', 'popupmenu', 'String', opts_str, ...
                    'Units', 'normalized', 'Position', [0.22 0.35 0.15 0.2], 'Value', 3, ...
                    'Callback', @updatePreview, 'FontSize', 11);

    % Dropdown for Atlas Dim 3
    uicontrol(controlPanel, 'Style', 'text', 'String', 'Map to Atlas Dim 3 (Z-axis / Slices):', ...
              'Units', 'normalized', 'Position', [0.02 0.1 0.2 0.2], ...
              'HorizontalAlignment', 'left', 'FontSize', 11, 'BackgroundColor', 'w');
    dd3 = uicontrol(controlPanel, 'Style', 'popupmenu', 'String', opts_str, ...
                    'Units', 'normalized', 'Position', [0.22 0.1 0.15 0.2], 'Value', 5, ...
                    'Callback', @updatePreview, 'FontSize', 11);

    % Status text (for warnings)
    statusText = uicontrol(controlPanel, 'Style', 'text', 'String', '', ...
                           'Units', 'normalized', 'Position', [0.4 0.35 0.3 0.3], ...
                           'ForegroundColor', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

    % Confirm Button
    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Confirm & Save Orientation', ...
              'Units', 'normalized', 'Position', [0.75 0.2 0.2 0.5], ...
              'FontSize', 12, 'FontWeight', 'bold', 'Callback', @confirmCallback);

    % --- 2. Top Images Panel (Using panel for clean grid layout) ---
    imagePanel = uipanel(fig, 'Position', [0.0 0.2 1.0 0.8], 'BorderType', 'none', 'BackgroundColor', 'w');
    
    p = panel(imagePanel);
    p.pack('v', 2); % 2 Rows (Top = Atlas, Bottom = Sample)
    p(1).pack('h', 5); % 5 Columns for Atlas slices
    p(2).pack('h', 5); % 5 Columns for Sample slices
    
    p.margin = [10 10 10 20]; % [left bottom right top]
    p.de.margin = 8;          % Margin between elements
    p.fontname = 'Arial';
    
    % Master Title using panel syntax
    p.title('\bfMatch Orientation: \rmTarget Atlas (Top) vs Sample (Bottom)');

    % Initialize Axes arrays
    ax_atlas = gobjects(1, 5);
    ax_sample = gobjects(1, 5);
    
    for i = 1:5
        % Select and format Atlas axes
        ax_atlas(i) = p(1, i).select();
        axis(ax_atlas(i), 'image', 'off');
        colormap(ax_atlas(i), gray);
        
        % Select and format Sample axes
        ax_sample(i) = p(2, i).select();
        axis(ax_sample(i), 'image', 'off');
        colormap(ax_sample(i), gray);
    end
    
    % Draw the fixed Atlas template
    Nz_atlas = size(atlasvol, 3);
    n_preview = min(5, Nz_atlas);
    atlas_preview_slices = round(linspace(Nz_atlas * 0.1, Nz_atlas * 0.9, n_preview));
    
    for i = 1:n_preview
        currim_atlas = atlasvol(:, :, atlas_preview_slices(i));
        % Improve contrast for visualization
        if max(currim_atlas(:)) > 0
            currim_atlas = adapthisteq(mat2gray(currim_atlas));
        end
        
        imshow(currim_atlas, [], 'Parent', ax_atlas(i));
        title(ax_atlas(i), sprintf('Atlas Z: %d', atlas_preview_slices(i)), 'FontSize', 11, 'Color', [0 0.4 0.7]);
        
        if i == 1
            ylabel(ax_atlas(i), 'TARGET ATLAS', 'FontWeight', 'bold', 'FontSize', 12, 'Visible', 'on', 'Color', [0 0.4 0.7]);
            set(ax_atlas(i), 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');
        end
    end

    % Initialize the preview for the sample
    updatePreview();

    % Wait for the user to click Confirm or close the window
    uiwait(fig);

    %% --- Nested Callback Functions ---
    
    function updatePreview(~, ~)
        % Get current selected values
        v1 = opts_val(dd1.Value);
        v2 = opts_val(dd2.Value);
        v3 = opts_val(dd3.Value);
        
        current_pvec = [v1, v2, v3];
        
        % Validate: dimensions must be unique (ignoring signs)
        if numel(unique(abs(current_pvec))) ~= 3
            statusText.String = 'INVALID: Dimensions must be unique!';
            for a = 1:5
                cla(ax_sample(a)); 
                title(ax_sample(a), 'Preview Disabled'); 
            end
            p.title('\bfMatch Orientation: \rm\color{red}Invalid Mapping Detected');
            return;
        else
            statusText.String = '';
            p.title({sprintf(...
                '\\bfMatch Orientation: \\rmTarget Atlas (Top) vs Sample (Bottom) | Current Mapping: [%d, %d, %d]', ...
                current_pvec(1), current_pvec(2), current_pvec(3)), ...
                ' '});
        end
        
        % Apply permutation
        tempvol = permuteBrainVolume(backvol, current_pvec);

        % Extract 5 slices along the new Dim 3
        Nz_new = size(tempvol, 3);
        sample_preview_slices = round(linspace(Nz_new * 0.1, Nz_new * 0.9, n_preview));
        
        % Update plots
        for idx = 1:n_preview
            islice = sample_preview_slices(idx);
            currim = tempvol(:, :, islice);
            
            % Enhance contrast dynamically
            currim = imadjust(currim);

            
            imshow(currim, [], 'Parent', ax_sample(idx));
            title(ax_sample(idx), sprintf('Sample Z: %d', islice), 'FontSize', 11, 'Color', [0.7 0.2 0]);
            
            if idx == 1
                ylabel(ax_sample(idx), 'YOUR SAMPLE', 'FontWeight', 'bold', 'FontSize', 12, 'Visible', 'on', 'Color', [0.7 0.2 0]);
                set(ax_sample(idx), 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');
            end
        end
    end

    function confirmCallback(~, ~)
        % Save the confirmed vector and resume execution
        v1 = opts_val(dd1.Value);
        v2 = opts_val(dd2.Value);
        v3 = opts_val(dd3.Value);
        
        current_pvec = [v1, v2, v3];
        
        if numel(unique(abs(current_pvec))) == 3
            permvec = current_pvec;
            uiresume(fig);
            delete(fig);
        else
            errordlg('Please fix the permutation mapping before confirming. Dimensions must be unique.', 'Invalid Mapping');
        end
    end

    function closeFigure(~, ~)
        % Handle user closing the window via the 'X' button
        selection = questdlg('Orientation not saved. Are you sure you want to exit?', ...
                             'Close Window', 'Yes', 'No', 'Yes');
        switch selection
            case 'Yes'
                uiresume(fig);
                delete(fig);
                error('Orientation mapping was cancelled by the user.');
            case 'No'
                return;
        end
    end
end