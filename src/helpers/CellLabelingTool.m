function CellLabelingTool()
% CELLLABELINGTOOL Interactive app to label cells (Numeric: 1=Cell, 0=Noise, NaN=Unsorted)
%
% Controls:
%   [C] : Label Cell (1)
%   [N] : Label Noise (0)
%   [U] : Unlabel (NaN)
%   [<-]: Back
%   [->]: Next / Skip
%   [S] : Save Progress

    % --- 1. File Selection ---
    [fName, fPath] = uigetfile('*.mat', 'Select Data File');
    if isequal(fName, 0), return; end
    fullPath = fullfile(fPath, fName);
    
    fprintf('Loading source data from %s...\n', fName);
    data = load(fullPath);
    
    if ~isfield(data, 'cell_images') 
        errordlg('File must contain "cell_images".', 'Data Error');
        return;
    end
    
    % Hardcoded sigma as requested
    data.sigmause = 6*[3 3 2]; 

    % --- 2. Preprocessing ---
    fprintf('Generating images... this may take a moment.\n');
    % Assumes prepareImagesForCNN is in your MATLAB path
    XTrain_Full = prepareImagesForCNN(data.cell_images, data.sigmause);
    N = size(XTrain_Full, 4);
    
    % --- 3. Session Setup ---
    [~, name, ~] = fileparts(fName);
    saveFile = fullfile(fPath, [name '_labeled.mat']);
    
    % Numeric Labels: NaN = Unsorted, 1 = Cell, 0 = Noise
    currentLabels = nan(N, 1);
    
    % Random permutation (seeded by filename)
    rng(sum(double(name))); 
    viewOrder = randperm(N);
    
    ptr = 1;

    % --- RESUME LOGIC ---
    if isfile(saveFile)
        savedData = load(saveFile, 'labels', 'labeledIndices', 'viewOrder');
        choice = questdlg('Found existing save file. Resume?', 'Resume', 'Yes', 'No', 'Yes');
        
        if strcmp(choice, 'Yes')
            if isfield(savedData, 'viewOrder'), viewOrder = savedData.viewOrder; end
            
            if isfield(savedData, 'labeledIndices') && isfield(savedData, 'labels')
                idxRestore = savedData.labeledIndices;
                lblRestore = savedData.labels;
                
                % Map saved values back to the full array
                valid = idxRestore <= N;
                currentLabels(idxRestore(valid)) = lblRestore(valid);
                indsdone = find(ismembc(viewOrder,savedData.labeledIndices));
                ptr      = indsdone(end);
            end
        end
    end


    % --- 4. GUI ---
    f = figure('Name', ['Labeling: ' fName], 'NumberTitle', 'off', ...
        'Color', 'w', 'KeyPressFcn', @keyHandler, ...
        'CloseRequestFcn', @onClose);
    
    ax = axes('Parent', f, 'Position', [0.05 0.15 0.9 0.75]);
    
    lblStatus = uicontrol(f, 'Style', 'text', 'FontSize', 14, ...
        'Units', 'normalized', 'Position', [0.05 0.02 0.9 0.08], ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'center');

    uicontrol(f, 'Style', 'text', 'FontSize', 10, 'ForegroundColor', [0.4 0.4 0.4], ...
        'Units', 'normalized', 'Position', [0.05 0.92 0.9 0.05], ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'center', ...
        'String', '[C] Cell (1)  |  [N] Noise (0)  |  [U] Unlabel (NaN)  |  [<-] Back  |  [->] Skip  |  [S] Save');

    updateDisplay();
    
    % --- Logic ---
    
    function updateDisplay()
        if ptr > N
            lblStatus.String = 'End of Dataset. Press S to Save.';
            cla(ax); return;
        end
        
        idx = viewOrder(ptr);
        imgBlock = XTrain_Full(:, :, :, idx);
        montageImg = [imgBlock(:,:,1), imgBlock(:,:,2), imgBlock(:,:,3)];
        
        % Calculate limits ignoring zeros (background)
        posVals = montageImg(montageImg > 0);
        if isempty(posVals)
            minVal = 0; maxVal = 1;
        else
            minVal = min(posVals, [], 'all');
            maxVal = max(montageImg, [], 'all');
        end
        
        imagesc(ax, montageImg, [minVal, maxVal]);
        colormap(ax, gray);       
        axis(ax, 'image');        
        axis(ax, 'off');          

        % Status Feedback
        val = currentLabels(idx);
        if isnan(val)
            statStr = "Unlabeled";
            col = 'k'; 
        elseif val == 1
            statStr = "CELL (1)";
            col = [0 0.7 0]; 
        elseif val == 0
            statStr = "NOISE (0)";
            col = [0.8 0 0]; 
        else
            statStr = "Unknown";
            col = 'k';
        end
        
        lblStatus.String = sprintf('Image %d / %d  [%s]', ptr, N, statStr);
        lblStatus.ForegroundColor = col;
        title(ax, sprintf('ID: %d', idx), 'Color', [0.5 0.5 0.5], 'FontSize', 10);
    end

    function keyHandler(~, ev)
        switch lower(ev.Key)
            case 'c', applyLabel(1);
            case 'n', applyLabel(0);
            case 'u', applyLabel(nan); % Doesn't advance automatically
            case 'leftarrow', move(-1);
            case 'rightarrow', move(1);
            case 's', saveData();
        end
    end

    function applyLabel(val)
        if ptr > N, return; end
        idx = viewOrder(ptr);
        currentLabels(idx) = val;
        
        if isnan(val)
            updateDisplay(); % Stay here if unlabeling
        else
            move(1); % Advance if labeling
        end
    end

    function move(dir)
        newPtr = ptr + dir;
        if newPtr >= 1 && newPtr <= N + 1
            ptr = newPtr;
            updateDisplay();
        end
    end

    function saveData(silent)
        if nargin < 1, silent = false; end
        
        % Filter: Select ONLY labeled entries (not NaN)
        labeledMask = ~isnan(currentLabels);
        
        XTrain = XTrain_Full(:, :, :, labeledMask);
        labels = currentLabels(labeledMask);
        labeledIndices = find(labeledMask);
        
        save(saveFile, 'XTrain', 'labels', 'labeledIndices', 'viewOrder', '-v7.3');
        
        if ~silent
            fprintf('Saved %d labeled images to %s\n', length(labels), saveFile);
            title(ax, '--- SAVED ---', 'Color', 'b');
            pause(0.2);
            updateDisplay();
        end
    end

    function onClose(~,~)
        selection = questdlg('Save changes?', 'Close', 'Yes', 'No', 'Cancel', 'Yes');
        switch selection
            case 'Yes', saveData(true); delete(f);
            case 'No', delete(f);
        end
    end
end