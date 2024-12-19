function volume_control_points_selector(vol1, vol2)
    % Create the UI figure
    fig = uifigure('Name', '3D Control Points Selector', 'Position', [100 100 1200 600]);

    % Data structure to store control points
    controlPoints = nan(100, 3, 2);

    % Grid layout for UI components
    gl = uigridlayout(fig, [1, 3], 'ColumnWidth', {'1x', '1x', 'fit'});

    % Panel for volume 1 viewer
    panel1 = uipanel(gl, 'Title', 'Volume 1');
    panel1.Layout.Column = 1;

    % Panel for volume 2 viewer
    panel2 = uipanel(gl, 'Title', 'Volume 2');
    panel2.Layout.Column = 2;

    % Initialize orthoslice viewers
    viewer1 = orthosliceViewer(vol1, 'Parent', panel1);
    viewer2 = orthosliceViewer(vol2, 'Parent', panel2);

    % Panel for controls
    setupControlPanel(gl);

    % Nested function to setup the control panel and buttons
    function setupControlPanel(gridLayout)
        controlPanel = uipanel(gridLayout, 'Title', 'Controls');
        controlPanel.Layout.Column = 3;
        glb = uigridlayout(controlPanel, [2, 1]);

        % Button to add point simultaneously from both viewers
        btnAddPoints = uibutton(glb, 'Text', 'Add Matched Points', 'ButtonPushedFcn', @(btn, event) addPoints());
        btnAddPoints.Layout.Row = 1;

        % Button to finish selection
        btnFinish = uibutton(glb, 'Text', 'Finish Selection', 'ButtonPushedFcn', @(btn, event) finishSelection());
        btnFinish.Layout.Row = 2;
    end

    % Function to add points from both orthosliceViewers simultaneously
    function addPoints()
        % Get current slice numbers and positions
        sliceNum1 = viewer1.SliceNumbers;
        sliceNum2 = viewer2.SliceNumbers;

        % Append to control points structure
        controlPoints(end+1, :, :) = cat(3, sliceNum1, sliceNum2);

        % Visual feedback: plot points on the slice viewers
        ax1 = viewer1.getAxesHandles;
        hold(ax1, 'on');
        plot3(ax1, sliceNum1(2), sliceNum1(1), sliceNum1(3), 'ro');
        hold(ax1, 'off');
    
        ax2 = viewer2.getAxesHandles;
        hold(ax2, 'on');
        plot3(ax2, sliceNum2(2), sliceNum2(1), sliceNum2(3), 'ro');
        hold(ax2, 'off');

        disp('Matched points added');
    end

    % Function to handle finishing the selection
    function finishSelection()
        disp('Selection finished. Control points:');
        disp('Volume 1 points:');
        disp(controlPoints(:,:,1));
        disp('Volume 2 points:');
        disp(controlPoints(:,:,2));
    end
end