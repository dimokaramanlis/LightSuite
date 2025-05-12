function InteractiveSliceReorder(optionalVolumePath)
% Interactively reorders slices from a single multi-page TIFF file.
% Uses a common _processing_decisions.txt file for order and flip states.
% Can take an optional full path to the TIFF volume as input.
% Implements a paged display for better viewing of many slices.
% Loads a previously saved order, visually reordering images in the GUI if specified.

    % --- Configurable Display Parameters ---
    displayRows = 2; % Number of rows of images per page
    displayCols = 3; % Number of columns of images per page
    itemsPerPage = displayRows * displayCols;

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

    originalSliceImagesData = cell(numSlices, 1); % Store of original images
    hWaitBar = waitbar(0, 'Loading slices from TIFF file...');
    try
        for i = 1:numSlices
            waitbar(i/numSlices, hWaitBar, sprintf('Loading slice %d/%d...', i, numSlices));
            originalSliceImagesData{i} = imread(inputFileFullPath, i, 'Info', tiffInfo);
        end
    catch ME
        if ishandle(hWaitBar); close(hWaitBar); end
        errordlg(['Error loading slices from TIFF: ' ME.message], 'Image Loading Error');
        return;
    end
    if ishandle(hWaitBar); close(hWaitBar); end

    [filePathStr, baseName, ~] = fileparts(inputFileFullPath);
    % Corrected filename construction to match the one used in SliceFlipper if it was changed there.
    % Using strcat as per user's note for consistency, though simple concatenation [baseName, '_suffix'] also works.
    processingDecisionsFilename = fullfile(filePathStr, strcat(baseName, '_processing_decisions.txt'));


    % --- GUI Creation ---
    screenSize = get(0, 'ScreenSize');
    guiWidthFraction = 0.8; 
    guiHeightFraction = 0.8;
    
    estTileSize = 150; 
    if numSlices < itemsPerPage && numSlices > 0
        if displayCols <=2 && displayRows <=2; estTileSize = 200; end
    end

    controlPanelRelHeight = 0.09; 

    figContentWidth = displayCols * estTileSize + displayCols * 10; 
    figContentHeight = displayRows * estTileSize + displayRows * 10;

    figWidth = min(screenSize(3)*guiWidthFraction, figContentWidth + 60); 
    figHeight = min(screenSize(4)*guiHeightFraction, figContentHeight / (1-controlPanelRelHeight) + 40);

    figPos = [(screenSize(3)-figWidth)/2, (screenSize(4)-figHeight)/2, figWidth, figHeight];

    guiFig = figure('Name', 'Interactive Slice Reorderer (Common Decisions)', ...
        'NumberTitle', 'off', 'Toolbar', 'none', 'Menubar', 'none', ...
        'Color', 'w', 'Units', 'pixels', 'Position', figPos, ...
        'CloseRequestFcn', @cancelAndClose); 

    mainPanel = uipanel('Parent', guiFig, 'Units', 'normalized', ...
                        'Position', [0 controlPanelRelHeight 1 (1-controlPanelRelHeight)], ...
                        'BorderType','none');

    controlPanel = uipanel('Parent', guiFig, 'Units', 'normalized', ...
                           'Position', [0 0 1 controlPanelRelHeight]); 

    guiData = struct();
    guiData.inputFileFullPath = inputFileFullPath; 
    guiData.originalSliceImagesData = originalSliceImagesData; 
    guiData.numSlices = numSlices;
    guiData.processingDecisionsFilename = processingDecisionsFilename;
    
    % Initialize state variables
    guiData.visualOrderOriginalIndices = (1:numSlices)'; 
    guiData.numberAssignedToOriginalSlice = nan(numSlices, 1); 
    guiData.flipState = zeros(numSlices, 1); % 0 for no flip, 1 for flip
    
    guiData.textHandles = gobjects(numSlices, 1); 
    
    guiData.displayRows = displayRows;
    guiData.displayCols = displayCols;
    guiData.itemsPerPage = itemsPerPage;
    guiData.currentPage = 1;
    guiData.totalPages = ceil(numSlices / itemsPerPage);
    if guiData.totalPages == 0; guiData.totalPages = 1; end 

    guiData.axesHandlesOnPage = gobjects(itemsPerPage, 1); 
    guiData.imageHandlesOnPage = gobjects(itemsPerPage, 1);
    guiData.currentTileOriginalIndices = nan(itemsPerPage, 1); 

    guiData.mainPanel = mainPanel;
    guiData.tileHandle = []; 
    guiData.estTileSize = estTileSize; 

    uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Cancel', ... 
        'Units', 'normalized', 'Position', [0.02 0.15 0.15 0.7], ... 
        'Callback', {@cancelAndCloseCallback, guiFig}, 'FontSize', 9);

    guiData.prevPageButton = uicontrol(controlPanel, 'Style', 'pushbutton', 'String', '< Prev', ...
        'Units', 'normalized', 'Position', [0.18 0.15 0.15 0.7], ... 
        'Callback', {@prevPageButtonCallback, guiFig}, 'FontSize', 9, 'Enable', 'off');
    
    guiData.pageInfoText = uicontrol(controlPanel, 'Style', 'text', 'String', sprintf('Page %d/%d', guiData.currentPage, guiData.totalPages), ...
        'Units', 'normalized', 'Position', [0.34 0.15 0.15 0.7], 'FontSize', 9, 'BackgroundColor', 'w');

    guiData.nextPageButton = uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Next >', ...
        'Units', 'normalized', 'Position', [0.50 0.15 0.15 0.7], ... 
        'Callback', {@nextPageButtonCallback, guiFig}, 'FontSize', 9);

    guiData.doneButton = uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Done & Save Decisions', ...
        'Units', 'normalized', 'Position', [0.68 0.15 0.30 0.7], ... 
        'Callback', {@doneButtonCallback, guiFig}, 'Enable', 'off', 'FontSize', 9); 
    
    guidata(guiFig, guiData); 

    % Attempt to load previously saved processing decisions
    if exist(guiData.processingDecisionsFilename, 'file') && numSlices > 0
        try
            % Corrected: Use readtable to handle headers
            loadedTable = readtable(guiData.processingDecisionsFilename, 'Delimiter', '\t', 'ReadVariableNames', true);
            
            if height(loadedTable) == numSlices && ...
               all(ismember({'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'}, loadedTable.Properties.VariableNames))
                
                % Load FlipState
                if isnumeric(loadedTable.FlipState)
                    guiData.flipState = loadedTable.FlipState;
                else % Handle if FlipState was read as cell/string (e.g. from older manual edits)
                    tempFlip = zeros(numSlices,1);
                    for k_fs=1:numSlices
                        if iscell(loadedTable.FlipState) && strcmpi(strtrim(string(loadedTable.FlipState{k_fs})), "1")
                            tempFlip(k_fs) = 1;
                        elseif isstring(loadedTable.FlipState) && strcmpi(strtrim(loadedTable.FlipState(k_fs)), "1")
                            tempFlip(k_fs) = 1;
                        end
                    end
                    guiData.flipState = tempFlip;
                    warning('FlipState column from %s was not purely numeric. Attempted conversion.', guiData.processingDecisionsFilename);
                end

                % Load NewOrderOriginalIndex (this column defines the visual order if complete)
                loadedNewOrderCol = loadedTable.NewOrderOriginalIndex; 
                
                % Populate numberAssignedToOriginalSlice
                tempNumAssigned = nan(numSlices,1);
                for k_new_pos = 1:numSlices % k_new_pos is the new order number
                    original_idx_at_k = loadedNewOrderCol(k_new_pos); % original slice index that should be at this new_pos
                    if ~isnan(original_idx_at_k) && original_idx_at_k >= 1 && original_idx_at_k <= numSlices
                        if isnan(tempNumAssigned(original_idx_at_k)) 
                            tempNumAssigned(original_idx_at_k) = k_new_pos;
                        else
                            warning('Duplicate original index %d found in loaded order file for new position. Using first occurrence.', original_idx_at_k);
                        end
                    end
                end
                guiData.numberAssignedToOriginalSlice = tempNumAssigned;

                % Determine visual order
                isLoadedOrderValidPermutation = ~any(isnan(loadedNewOrderCol)) && ...
                                                length(unique(loadedNewOrderCol)) == numSlices && ...
                                                all(ismember(loadedNewOrderCol, 1:numSlices));
                if isLoadedOrderValidPermutation
                    guiData.visualOrderOriginalIndices = loadedNewOrderCol;
                else
                    guiData.visualOrderOriginalIndices = (1:numSlices)'; 
                    if any(~isnan(loadedNewOrderCol)) % Only warn if there was some attempt at ordering
                        disp('Loaded order column was not a full permutation. Displaying slices in original sequence initially.');
                    end
                end
                
                guidata(guiFig, guiData); 
                disp(['Loaded decisions from: ', guiData.processingDecisionsFilename]);
            else
                warning('Dimension or column name mismatch in %s. Expected %dx3 table with specific headers. Starting fresh.', guiData.processingDecisionsFilename, numSlices);
            end
        catch ME_load
            warning('Error loading/parsing %s: %s. Starting fresh.', guiData.processingDecisionsFilename, ME_load.message);
        end
    end
    
    displayCurrentPage(guiFig); 
    updateDoneButtonState(guiFig); 
    updatePageButtons(guiFig);
end

% --- Page Display Function ---
function displayCurrentPage(fig) 
    guiData = guidata(fig);
    estTileSize = guiData.estTileSize; 

    if ~isempty(guiData.tileHandle) && ishandle(guiData.tileHandle)
        delete(guiData.tileHandle); 
    else 
        childHandles = allchild(guiData.mainPanel); 
        if ~isempty(childHandles); delete(childHandles); end
    end
    
    guiData.tileHandle = tiledlayout(guiData.mainPanel, guiData.displayRows, guiData.displayCols, ...
                                     'TileSpacing', 'compact', 'Padding', 'compact'); 

    guiData.axesHandlesOnPage = gobjects(guiData.itemsPerPage, 1);
    guiData.imageHandlesOnPage = gobjects(guiData.itemsPerPage, 1);
    guiData.currentTileOriginalIndices(:) = NaN; 

    startIndexOverallVisual = (guiData.currentPage - 1) * guiData.itemsPerPage + 1;
    
    tileCounter = 0;
    for i_visual_slot_on_page = 1:guiData.itemsPerPage
        tileCounter = tileCounter + 1;
        currentOverallVisualPosition = startIndexOverallVisual + i_visual_slot_on_page - 1;

        if currentOverallVisualPosition > guiData.numSlices; break; end 

        originalSliceIdxToDisplay = guiData.visualOrderOriginalIndices(currentOverallVisualPosition);
        if isnan(originalSliceIdxToDisplay) 
            ax = nexttile(guiData.tileHandle); 
            guiData.axesHandlesOnPage(tileCounter) = ax;
            text(ax, 0.5, 0.5, 'Empty Slot', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            guiData.imageHandlesOnPage(tileCounter) = gobjects(1);
            guiData.currentTileOriginalIndices(tileCounter) = NaN; 
            continue;
        end
        
        guiData.currentTileOriginalIndices(tileCounter) = originalSliceIdxToDisplay;

        ax = nexttile(guiData.tileHandle);
        guiData.axesHandlesOnPage(tileCounter) = ax;
        
        currentImage = guiData.originalSliceImagesData{originalSliceIdxToDisplay};
        if guiData.flipState(originalSliceIdxToDisplay) == 1
            currentImage = fliplr(currentImage);
        end

        if isempty(currentImage) 
            text(ax, 0.5, 0.5, sprintf('Orig. Slice %d Empty', originalSliceIdxToDisplay), 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            guiData.imageHandlesOnPage(tileCounter) = gobjects(1); 
            continue;
        end

        [imgH, imgW, ~] = size(currentImage);
        avgDim = mean([imgH, imgW]);
        if avgDim == 0; avgDim = estTileSize; end
        resizeFactor = estTileSize / avgDim;
        if resizeFactor > 1; resizeFactor = 1; end 
        if resizeFactor <=0; resizeFactor = 0.1; end 

        if ismatrix(currentImage)
            imgHandle = image(ax, imresize(currentImage, resizeFactor));
            colormap(ax, gray);
        else
            imgHandle = image(ax, imresize(currentImage, resizeFactor));
        end
        axis(ax, 'image', 'off');
        set(imgHandle, 'ButtonDownFcn', {@clickSliceCallback, fig, tileCounter}); 
        guiData.imageHandlesOnPage(tileCounter) = imgHandle;

        newOrderNumberToShow = guiData.numberAssignedToOriginalSlice(originalSliceIdxToDisplay);
        if ~isnan(newOrderNumberToShow)
            if ishandle(guiData.textHandles(newOrderNumberToShow))
                delete(guiData.textHandles(newOrderNumberToShow));
            end
            guiData.textHandles(newOrderNumberToShow) = text(ax, 'Units', 'normalized', ...
                'Position', [0.05 0.15], 'String', num2str(newOrderNumberToShow), ...
                'FontSize', 16, 'Color', 'white', 'FontWeight', 'bold', ...
                'BackgroundColor', [0 0 0 0.5]);
        end
    end
    
    numAssigned = sum(~isnan(guiData.numberAssignedToOriginalSlice));
    titleString = sprintf('Click to assign order. %d of %d slices ordered.', numAssigned, guiData.numSlices);
    title(guiData.tileHandle, titleString, 'FontSize', 12, 'FontWeight', 'bold');

    guidata(fig, guiData);
end


% --- Callbacks ---
function clickSliceCallback(~, ~, fig, tileNumOnPage)
    guiData = guidata(fig);
    originalSliceIdxClicked = guiData.currentTileOriginalIndices(tileNumOnPage); 

    if isnan(originalSliceIdxClicked)
        disp('Clicked on an empty or unmapped tile.');
        return;
    end

    currentAssignedNewOrder = guiData.numberAssignedToOriginalSlice(originalSliceIdxClicked);

    if ~isnan(currentAssignedNewOrder) 
        if ishandle(guiData.textHandles(currentAssignedNewOrder))
            delete(guiData.textHandles(currentAssignedNewOrder));
            guiData.textHandles(currentAssignedNewOrder) = gobjects(1); 
        end
        guiData.numberAssignedToOriginalSlice(originalSliceIdxClicked) = NaN;
    else 
        nextAvailableNewOrderNumber = NaN;
        for k_new_order = 1:guiData.numSlices
            isNumUsed = false;
            for orig_idx_check = 1:guiData.numSlices
                if guiData.numberAssignedToOriginalSlice(orig_idx_check) == k_new_order
                    isNumUsed = true;
                    break;
                end
            end
            if ~isNumUsed
                nextAvailableNewOrderNumber = k_new_order;
                break;
            end
        end
        
        if isnan(nextAvailableNewOrderNumber)
            disp('All order numbers seem assigned. Cannot assign new one.');
            return; 
        end
        
        guiData.numberAssignedToOriginalSlice(originalSliceIdxClicked) = nextAvailableNewOrderNumber;
        axForText = guiData.axesHandlesOnPage(tileNumOnPage);
        if ishandle(guiData.textHandles(nextAvailableNewOrderNumber))
            delete(guiData.textHandles(nextAvailableNewOrderNumber));
        end
        guiData.textHandles(nextAvailableNewOrderNumber) = text(axForText, 'Units', 'normalized', ...
            'Position', [0.05 0.15], 'String', num2str(nextAvailableNewOrderNumber), ...
            'FontSize', 16, 'Color', 'white', 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.5]); 
    end
    
    guidata(fig, guiData); 
    displayCurrentPage(fig); 
    updateDoneButtonState(fig); 
end

function prevPageButtonCallback(~, ~, fig) 
    guiData = guidata(fig);
    if guiData.currentPage > 1
        guiData.currentPage = guiData.currentPage - 1;
        guidata(fig, guiData);
        displayCurrentPage(fig); 
        updatePageButtons(fig);
    end
end

function nextPageButtonCallback(~, ~, fig) 
    guiData = guidata(fig);
    if guiData.currentPage < guiData.totalPages
        guiData.currentPage = guiData.currentPage + 1;
        guidata(fig, guiData);
        displayCurrentPage(fig); 
        updatePageButtons(fig);
    end
end

function doneButtonCallback(~, ~, fig) 
    guiData = guidata(fig);
    all_numbers_assigned = true;
    if guiData.numSlices > 0
        for k_num = 1:guiData.numSlices
            if ~any(guiData.numberAssignedToOriginalSlice == k_num)
                all_numbers_assigned = false;
                break;
            end
        end
    end

    if ~all_numbers_assigned && guiData.numSlices > 0
        warndlg('Not all slices have been assigned an order number from 1 to N.', 'Ordering Incomplete');
        return;
    end
    
    choice = questdlg('Are you sure you want to save the processing decisions?', ...
        'Confirm Save', 'Yes, Save', 'No', 'No');
    
    if strcmp(choice, 'Yes, Save')
        try
            saveProcessingDecisions(guiData); 
            msgbox('Processing decisions have been successfully saved.', 'Success', 'modal');
        catch ME
            errordlg(['Error saving decisions file: ' ME.message], 'Saving Error');
        end
        delete(fig); 
    end
end

function cancelAndClose(fig, ~) 
    choice = questdlg('Are you sure you want to cancel? Changes will be lost if not saved.', ...
        'Confirm Cancel', 'Yes, Cancel', 'No', 'No');
    if strcmp(choice, 'Yes, Cancel')
        delete(fig);
    end
end

function cancelAndCloseCallback(~, ~, fig)
    cancelAndClose(fig); 
end

% --- Utility Functions ---
function updateTitle(fig) 
    guiData = guidata(fig);
    if ~isempty(guiData.tileHandle) && ishandle(guiData.tileHandle)
        numAssigned = sum(~isnan(guiData.numberAssignedToOriginalSlice)); 
        titleString = sprintf('Click to assign order. %d of %d slices ordered.', numAssigned, guiData.numSlices);
        title(guiData.tileHandle, titleString, 'FontSize', 12, 'FontWeight', 'bold');
    end
end

function updateDoneButtonState(fig)
    guiData = guidata(fig);
    all_numbers_assigned = true;
    if guiData.numSlices > 0
        for k_num = 1:guiData.numSlices
            if ~any(guiData.numberAssignedToOriginalSlice == k_num) 
                all_numbers_assigned = false;
                break;
            end
        end
    end

    if (all_numbers_assigned && guiData.numSlices > 0) || (guiData.numSlices == 0)
        set(guiData.doneButton, 'Enable', 'on');
    else
        set(guiData.doneButton, 'Enable', 'off');
    end
end

function updatePageButtons(fig)
    guiData = guidata(fig);
    set(guiData.pageInfoText, 'String', sprintf('Page %d/%d', guiData.currentPage, guiData.totalPages));
    if guiData.currentPage == 1
        set(guiData.prevPageButton, 'Enable', 'off');
    else
        set(guiData.prevPageButton, 'Enable', 'on');
    end
    if guiData.currentPage == guiData.totalPages || guiData.numSlices == 0
        set(guiData.nextPageButton, 'Enable', 'off');
    else
        set(guiData.nextPageButton, 'Enable', 'on');
    end
end

function saveProcessingDecisions(guiData)
    decisionsMatrix = zeros(guiData.numSlices, 3);
    if guiData.numSlices > 0
        decisionsMatrix(:,1) = (1:guiData.numSlices)'; 
        decisionsMatrix(:,2) = guiData.flipState;      
    
        newOrderColumn = nan(guiData.numSlices, 1);
        for original_idx = 1:guiData.numSlices
            assigned_new_order_num = guiData.numberAssignedToOriginalSlice(original_idx);
            if ~isnan(assigned_new_order_num)
                if assigned_new_order_num >= 1 && assigned_new_order_num <= guiData.numSlices
                     if isnan(newOrderColumn(assigned_new_order_num))
                        newOrderColumn(assigned_new_order_num) = original_idx;
                     else
                        warning('Attempting to assign original slice %d to new position %d, but it is already occupied by original slice %d. Check logic.', ...
                                original_idx, assigned_new_order_num, newOrderColumn(assigned_new_order_num));
                     end
                else
                    warning('Invalid assigned new order number %d for original slice %d. Skipping.', assigned_new_order_num, original_idx);
                end
            end
        end
        decisionsMatrix(:,3) = newOrderColumn;
    else % Case for 0 slices
        decisionsMatrix = zeros(0,3); % Ensure it's an empty matrix with 3 columns
    end
    
    hWaitBarSave = waitbar(0, 'Saving processing decisions...');
    try
        if guiData.numSlices > 0
            T = array2table(decisionsMatrix, 'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
            writetable(T, guiData.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t');
        else
            T_empty = table('Size', [0,3], 'VariableTypes', {'double', 'double', 'double'}, ...
                            'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
            writetable(T_empty, guiData.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t');
        end
        
        if ishandle(hWaitBarSave); close(hWaitBarSave); end
        disp(['Processing decisions saved to: ', guiData.processingDecisionsFilename]);
    catch ME
        if ishandle(hWaitBarSave); close(hWaitBarSave); end
        rethrow(ME); 
    end
end
