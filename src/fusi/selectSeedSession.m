function iseed = selectSeedSession(volumes, savepath, sessionnames, voxelsize_mm)
%SELECTSEEDSESSION Pick the session that will define the common anatomy space.
%   iseed = selectSeedSession(volumes) opens a figure where every session is
%   shown as a column and a few coronal slices (the last dimension of each
%   volume) are stacked as rows within that column. Clicking anywhere in a
%   column selects that session; 'Confirm' returns its index.
%
%   iseed = selectSeedSession(volumes, savepath) additionally caches the
%   choice in <savepath>/seed_session_for_anatomy.txt. If that file already
%   exists the GUI is skipped and the stored index is returned, so the whole
%   anatomy pipeline can be re-run without asking again.
%
%   Inputs:
%     volumes        - 4-D array [d1 d2 d3 Nsessions] of per-session volumes
%                      (typically the median across time of each recording).
%                      Coronal slices run along dim 3.
%     savepath       - (optional) folder holding seed_session_for_anatomy.txt.
%                      Pass '' to disable caching and always show the GUI.
%     sessionnames   - (optional) cellstr of Nsessions labels for the columns.
%     voxelsize_mm   - (optional) [d1 d2 d3] voxel size of the volumes, used
%                      only to draw each coronal slice with the correct aspect
%                      ratio. Defaults to isotropic.
%
%   Output:
%     iseed          - index into the 4th dimension of volumes.

Nsessions = size(volumes, 4);
Nz        = size(volumes, 3);

if nargin < 2, savepath = ''; end
if nargin < 3 || isempty(sessionnames)
    sessionnames = arrayfun(@(x) sprintf('session %d', x), 1:Nsessions, 'UniformOutput', false);
end
if nargin < 4 || isempty(voxelsize_mm)
    voxelsize_mm = [1 1 1];
end
voxelsize_mm = double(voxelsize_mm(:)).';
if isscalar(voxelsize_mm), voxelsize_mm = voxelsize_mm([1 1 1]); end
%==========================================================================
% a previous choice short-circuits everything
seedfile = '';
if ~isempty(savepath)
    seedfile = fullfile(savepath, 'seed_session_for_anatomy.txt');
    fprintf('Looking for seed session in %s\n', seedfile)
    if exist(seedfile, 'file')
        iseed = readmatrix(seedfile);
        iseed = iseed(1);
        if iseed < 1 || iseed > Nsessions || mod(iseed, 1) ~= 0
            error('selectSeedSession:badCache', ...
                'Stored seed session %g is not a valid index into %d sessions (%s).', ...
                iseed, Nsessions, seedfile);
        end
        fprintf('Found it, seed session is %d (%s)\n', iseed, sessionnames{iseed})
        return;
    end
    fprintf('Not found, you have to pick a seed session, check GUI\n')
end
%==========================================================================
if Nsessions == 1
    iseed = 1;
else
    iseed = runSelectionGui(volumes, sessionnames, voxelsize_mm, Nz);
end
%==========================================================================
if ~isempty(seedfile)
    writematrix(iseed, seedfile);
    fprintf('Seed session saved as %d to %s\n', iseed, seedfile);
end
%==========================================================================
end

function iseed = runSelectionGui(volumes, sessionnames, voxelsize_mm, Nz)
%RUNSELECTIONGUI Sessions as columns, coronal slices as rows, click to pick.

Nsessions = numel(sessionnames);
Nshow     = min(6, Nz);
islices   = unique(round(linspace(0.15*Nz, 0.85*Nz, Nshow)));
islices   = max(min(islices, Nz), 1);
Nshow     = numel(islices);

% Coronal slices are [dim1 x dim2] images. DataAspectRatio uses inverse
% spacings: screen-length-per-unit is proportional to 1/DataAspectRatio.
slice_aspect = [1/voxelsize_mm(2), 1/voxelsize_mm(1), 1];

% One contrast range per session so slices within a column are comparable.
climits = nan(Nsessions, 2);
for isess = 1:Nsessions
    currvol = volumes(:, :, islices, isess);
    climits(isess, :) = quantile(single(currvol(:)), [0.05 0.95]);
    if ~(climits(isess, 2) > climits(isess, 1))
        climits(isess, :) = [0 1];
    end
end

selected  = 1;
confirmed = false;

fig = figure('Name', 'Select Seed Session for Anatomy', ...
             'Position', [80 80 min(220*Nsessions + 120, 1800) 820], ...
             'MenuBar', 'none', 'ToolBar', 'none', ...
             'Color', 'w', 'NumberTitle', 'off', ...
             'CloseRequestFcn', @closeFigure, 'KeyPressFcn', @keyCallback);

% --- 1. Bottom control panel -------------------------------------------
controlPanel = uipanel(fig, 'Position', [0.0 0.0 1.0 0.11], ...
                       'BackgroundColor', 'w', 'BorderType', 'none');

uicontrol(controlPanel, 'Style', 'text', 'Units', 'normalized', ...
          'Position', [0.02 0.55 0.5 0.35], 'HorizontalAlignment', 'left', ...
          'FontSize', 11, 'BackgroundColor', 'w', ...
          'String', 'Click a column to select that session (arrow keys also work, Enter confirms).');

statusText = uicontrol(controlPanel, 'Style', 'text', 'Units', 'normalized', ...
                       'Position', [0.02 0.12 0.5 0.35], 'HorizontalAlignment', 'left', ...
                       'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', ...
                       'ForegroundColor', [0.7 0.2 0]);

uicontrol(controlPanel, 'Style', 'pushbutton', 'String', 'Confirm & Save Seed Session', ...
          'Units', 'normalized', 'Position', [0.75 0.2 0.22 0.6], ...
          'FontSize', 12, 'FontWeight', 'bold', 'Callback', @confirmCallback);

% --- 2. One panel per session, filled with coronal slices ---------------
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
          'Position', [0.0 0.955 1.0 0.035], 'BackgroundColor', 'w', ...
          'FontSize', 13, 'FontWeight', 'bold', ...
          'String', 'Seed session: pick the column with the cleanest, most complete brain coverage');

colpanel = gobjects(1, Nsessions);
xmargin  = 0.005;
colwidth = (1 - 2*xmargin)/Nsessions;
for isess = 1:Nsessions
    colpanel(isess) = uipanel(fig, ...
        'Position', [xmargin + (isess-1)*colwidth, 0.12, colwidth, 0.83], ...
        'Title', sprintf('%d: %s', isess, sessionnames{isess}), ...
        'TitlePosition', 'centertop', 'FontSize', 10, ...
        'BackgroundColor', 'w', 'BorderType', 'line', ...
        'ButtonDownFcn', makeSelector(isess));

    for islice = 1:Nshow
        ax = axes('Parent', colpanel(isess), 'Units', 'normalized', ...
            'Position', [0.02, 1 - islice/Nshow + 0.01, 0.96, 1/Nshow - 0.02]);
        imh = imagesc(ax, volumes(:, :, islices(islice), isess));
        colormap(ax, gray);
        axis(ax, 'image');
        set(ax, 'CLim', climits(isess, :), 'DataAspectRatio', slice_aspect, ...
            'XTick', [], 'YTick', [], 'Box', 'off', ...
            'XColor', 'none', 'YColor', 'none', ...
            'ButtonDownFcn', makeSelector(isess));
        set(imh, 'ButtonDownFcn', makeSelector(isess));

        if isess == 1
            text(ax, -0.04, 0.5, sprintf('z = %d', islices(islice)), ...
                'Units', 'normalized', 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 9, 'Color', [0.3 0.3 0.3]);
        end
    end
end

updateHighlight();
uiwait(fig);

if ~confirmed
    error('selectSeedSession:cancelled', 'Seed session selection was cancelled by the user.');
end
iseed = selected;

%% --- Nested callbacks -------------------------------------------------

    function fcn = makeSelector(isess)
        fcn = @(~, ~) selectColumn(isess);
    end

    function selectColumn(isess)
        selected = isess;
        updateHighlight();
        % double-click on a column is a shortcut for confirming
        if strcmp(get(fig, 'SelectionType'), 'open')
            confirmCallback();
        end
    end

    function updateHighlight(~, ~)
        for ii = 1:Nsessions
            if ii == selected
                set(colpanel(ii), 'HighlightColor', [0.7 0.2 0], 'BorderWidth', 4, ...
                    'ForegroundColor', [0.7 0.2 0], 'FontWeight', 'bold');
            else
                set(colpanel(ii), 'HighlightColor', [0.8 0.8 0.8], 'BorderWidth', 1, ...
                    'ForegroundColor', [0.3 0.3 0.3], 'FontWeight', 'normal');
            end
        end
        statusText.String = sprintf('Selected: %d (%s)', selected, sessionnames{selected});
    end

    function keyCallback(~, evt)
        switch evt.Key
            case {'rightarrow', 'downarrow'}
                selected = min(selected + 1, Nsessions);
                updateHighlight();
            case {'leftarrow', 'uparrow'}
                selected = max(selected - 1, 1);
                updateHighlight();
            case {'return', 'space'}
                confirmCallback();
        end
    end

    function confirmCallback(~, ~)
        confirmed = true;
        uiresume(fig);
        delete(fig);
    end

    function closeFigure(~, ~)
        selection = questdlg('No seed session saved. Are you sure you want to exit?', ...
                             'Close Window', 'Yes', 'No', 'Yes');
        if strcmp(selection, 'Yes')
            uiresume(fig);
            delete(fig);
        end
    end
end
