function annotateGRINLens(savepath, varargin)
% ANNOTATEGRINLENS  GUI for annotating an optical fiber / GRIN lens position
%   in a LightSuite registered volume.
%
%   Loads the 20 um registration volume, lets the user navigate through slices
%   and click on the edges of the fiber.  Pressing F transforms the marked
%   points to Allen-CCF atlas space, fits a circle that represents the bottom
%   face of the cylinder, and opens a figure with atlas-region colour maps at
%   0, 50, 100, 150, 200 and 250 um from the lens bottom.
%
%   Usage:
%       annotateGRINLens(savepath)
%       annotateGRINLens(savepath, 'Diameter', 500)
%       annotateGRINLens(savepath, 'Diameter', 500, 'Channel', 1)
%
%   Required input:
%       savepath   – LightSuite folder that already contains
%                   regopts.mat and transform_params.mat
%
%   Optional name-value inputs:
%       Diameter   – fiber diameter in µm  (default 500)
%       Channel    – which registration channel to display (default 1)
%
%   Keyboard controls:
%       ←  / →       navigate slices
%       1 / 2 / 3    change slice dimension (AP / DV / ML)
%       Click        add edge point on current slice
%       Backspace    delete last added point
%       C            clear all points from current slice
%       S            save points to grin_fiber_points.mat
%       F            fit circle + save + open atlas images
%       Enter        jump to a specific slice number

p = inputParser;
addRequired(p, 'savepath');
addParameter(p, 'Diameter', 500, @isnumeric);
addParameter(p, 'Channel',  1,   @isnumeric);
parse(p, savepath, varargin{:});

%----------------------------------------------------------------------
% 1. Load registration data
%----------------------------------------------------------------------
optsfile = fullfile(savepath, 'regopts.mat');
trfile   = fullfile(savepath, 'transform_params.mat');
if ~exist(optsfile, 'file') || ~exist(trfile, 'file')
    error('annotateGRINLens: regopts.mat or transform_params.mat not found in %s', savepath);
end

opts_data = load(optsfile);
regopts   = opts_data.opts;
trstruct  = load(trfile);

% Locate registration volume (one TIFF per channel)
regvolpath = regopts.regvolpath;
if iscell(regvolpath)
    regvolpath = regvolpath{p.Results.Channel};
end
if ~exist(regvolpath, 'file')
    error('annotateGRINLens: registration volume not found at %s', regvolpath);
end

fprintf('Loading registration volume from %s ...\n', regvolpath);
volraw = readDownStack(regvolpath);

% Permute to atlas orientation (same convention as registration pipeline)
how_to_perm = trstruct.how_to_perm;
volperm     = permuteBrainVolume(volraw, how_to_perm);

% Normalise to uint8 for display
irand  = randperm(numel(volperm), min(numel(volperm), 1e5));
factv  = 255 / single(quantile(single(volperm(irand)), 0.999));
voldisp = uint8(min(single(volperm) * factv, 255));

fprintf('Volume size (permuted): %d x %d x %d\n', size(voldisp,1), size(voldisp,2), size(voldisp,3));

%----------------------------------------------------------------------
% 2. Initialise GUI state
%----------------------------------------------------------------------
gui_data = struct;
gui_data.savepath      = savepath;
gui_data.trstruct      = trstruct;
gui_data.regopts       = regopts;
gui_data.vol           = voldisp;
gui_data.how_to_perm   = how_to_perm;
gui_data.diameter_um   = p.Results.Diameter;
gui_data.registres     = regopts.registres;
gui_data.ori_pxsize    = trstruct.ori_pxsize;
gui_data.slice_dim     = 1;                             % navigate along dim 1 by default
gui_data.curr_slice    = round(size(voldisp, 1) / 2);
gui_data.points        = zeros(0, 3);                   % Nx3 [r c s] in permuted vol
gui_data.save_file     = fullfile(savepath, 'grin_fiber_points.mat');

% Re-load previously saved points if present
if exist(gui_data.save_file, 'file')
    sv = load(gui_data.save_file, 'points', 'diameter_um', 'slice_dim');
    if isfield(sv, 'points') && ~isempty(sv.points)
        gui_data.points = sv.points;
        fprintf('Loaded %d saved points.\n', size(sv.points, 1));
    end
    if isfield(sv, 'diameter_um');  gui_data.diameter_um = sv.diameter_um; end
    if isfield(sv, 'slice_dim');    gui_data.slice_dim   = sv.slice_dim;   end
    gui_data.curr_slice = round(size(voldisp, gui_data.slice_dim) / 2);
end

%----------------------------------------------------------------------
% 3. Build figure and layout
%----------------------------------------------------------------------
screen_sz = get(0, 'screensize');
fig_w = screen_sz(3) * 0.55;
fig_h = fig_w * 0.7;
gui_fig = figure( ...
    'Name',             'GRIN Lens Annotation', ...
    'KeyPressFcn',      @keypress, ...
    'WindowScrollWheelFcn', @scroll_slice, ...
    'Toolbar',          'none', ...
    'Menubar',          'none', ...
    'Color',            'w', ...
    'Units',            'pixels', ...
    'Position',         [(screen_sz(3)-fig_w)/2, (screen_sz(4)-fig_h)/2, fig_w, fig_h], ...
    'CloseRequestFcn',  @close_gui);

gui_data.pp = panel();
gui_data.pp.pack('v', {0.92, 0.08});
gui_data.pp.margin = [2 2 2 28];

gui_data.base_title = { ...
    ['\bfGRIN Lens Annotation  (diameter: ' num2str(gui_data.diameter_um) ' µm)'], ...
    ['\bfControls:\rm  ←/→ or scroll: slice  |  Click: add point  |  ' ...
     'Backspace: delete  |  C: clear slice  |  1/2/3: dim  |  F: fit atlas  |  S: save']};
gui_data.pp.title(gui_data.base_title);
gui_data.pp.fontname = 'Arial';

% Main image axis
ax = gui_data.pp(1).select();
gui_data.ax = ax;
ax.YDir     = 'reverse';
ax.Colormap = gray;
hold(ax, 'on');
axis(ax, 'image', 'off');

curr_im       = get_slice(gui_data);
gui_data.im_h = imagesc(curr_im, 'Parent', ax, 'ButtonDownFcn', @mouseclick);
clim(ax, [0 255]);

gui_data.pts_h = plot(ax, nan, nan, '.r', 'MarkerSize', 16, 'PickableParts', 'none');

% Status line
stat_ax = gui_data.pp(2).select();
axis(stat_ax, 'off');
gui_data.status_txt = text(stat_ax, 0.01, 0.5, '', ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'middle');

guidata(gui_fig, gui_data);
update_display(gui_fig);

end  % annotateGRINLens

%==========================================================================
% CALLBACKS
%==========================================================================

function keypress(gui_fig, event)
    gui_data    = guidata(gui_fig);
    needs_update = true;

    switch event.Key
        case 'leftarrow'
            gui_data.curr_slice = max(1, gui_data.curr_slice - 1);

        case 'rightarrow'
            max_sl = size(gui_data.vol, gui_data.slice_dim);
            gui_data.curr_slice = min(max_sl, gui_data.curr_slice + 1);

        case {'1', 'numpad1'}
            gui_data.slice_dim  = 1;
            gui_data.curr_slice = round(size(gui_data.vol, 1) / 2);

        case {'2', 'numpad2'}
            gui_data.slice_dim  = 2;
            gui_data.curr_slice = round(size(gui_data.vol, 2) / 2);

        case {'3', 'numpad3'}
            gui_data.slice_dim  = 3;
            gui_data.curr_slice = round(size(gui_data.vol, 3) / 2);

        case 'return'
            max_sl = size(gui_data.vol, gui_data.slice_dim);
            ans_   = inputdlg(sprintf('Jump to slice (1 – %d):', max_sl), 'Go to slice');
            if ~isempty(ans_)
                v = str2double(ans_{1});
                if ~isnan(v) && v >= 1 && v <= max_sl
                    gui_data.curr_slice = round(v);
                end
            end

        case 'backspace'
            if ~isempty(gui_data.points)
                gui_data.points(end, :) = [];
            end

        case 'c'
            idim = gui_data.slice_dim;
            keep = gui_data.points(:, idim) ~= gui_data.curr_slice;
            gui_data.points = gui_data.points(keep, :);

        case 's'
            guidata(gui_fig, gui_data);
            do_save(gui_data);
            flash_title(gui_data, '\bfSaved!');
            needs_update = false;

        case 'f'
            guidata(gui_fig, gui_data);
            do_save(gui_data);
            fit_and_show_atlas(gui_fig);
            return;

        otherwise
            needs_update = false;
    end

    if needs_update
        guidata(gui_fig, gui_data);
        update_display(gui_fig);
    end
end

function scroll_slice(gui_fig, event)
    gui_data = guidata(gui_fig);
    step     = event.VerticalScrollCount;
    max_sl   = size(gui_data.vol, gui_data.slice_dim);
    gui_data.curr_slice = max(1, min(max_sl, gui_data.curr_slice + step));
    guidata(gui_fig, gui_data);
    update_display(gui_fig);
end

function mouseclick(src, event)
    if ~strcmp(get(ancestor(src,'figure'), 'SelectionType'), 'normal')
        return;
    end
    gui_data = guidata(src);

    idim      = gui_data.slice_dim;
    alldims   = 1:3;
    free_dims = alldims(alldims ~= idim);  % the two image dimensions
    pt3d      = zeros(1, 3);

    % IntersectionPoint is [x y z] = [col row depth] in axes data coords
    xy = event.IntersectionPoint(1:2);        % [col, row]
    pt3d(free_dims(1)) = xy(2);               % row → first free dim
    pt3d(free_dims(2)) = xy(1);               % col → second free dim
    pt3d(idim)         = gui_data.curr_slice;

    gui_data.points = [gui_data.points; pt3d];
    guidata(src, gui_data);
    update_display(ancestor(src, 'figure'));
end

function close_gui(gui_fig, ~)
    gui_data = guidata(gui_fig);
    if ~isempty(gui_data.points)
        ch = questdlg('Save points before closing?', 'Close', 'Yes', 'No', 'Yes');
        if strcmp(ch, 'Yes')
            do_save(gui_data);
        end
    end
    delete(gui_fig);
end

%==========================================================================
% DISPLAY
%==========================================================================

function update_display(gui_fig)
    gui_data = guidata(gui_fig);

    % Update slice image
    curr_im = get_slice(gui_data);
    try; curr_im = adapthisteq(curr_im); catch; end
    set(gui_data.im_h, 'CData', curr_im);

    % Rescale axis to new image size (handles dim switches)
    ax = gui_data.ax;
    set(ax, 'XLim', [0.5, size(curr_im,2)+0.5], 'YLim', [0.5, size(curr_im,1)+0.5]);

    % Points on current slice
    idim      = gui_data.slice_dim;
    alldims   = 1:3;
    free_dims = alldims(alldims ~= idim);

    on_slice = gui_data.points(:, idim) == gui_data.curr_slice;
    pts_vis  = gui_data.points(on_slice, :);

    if ~isempty(pts_vis)
        set(gui_data.pts_h, ...
            'XData', pts_vis(:, free_dims(2)), ...
            'YData', pts_vis(:, free_dims(1)));
    else
        set(gui_data.pts_h, 'XData', nan, 'YData', nan);
    end

    % Status
    max_sl = size(gui_data.vol, idim);
    npts   = size(gui_data.points, 1);
    dim_labels = {'AP','DV','ML'};
    msg = sprintf('Slice %d / %d  (%s, dim %d)  |  %d total points  |  Diameter: %d µm', ...
        gui_data.curr_slice, max_sl, dim_labels{idim}, idim, npts, gui_data.diameter_um);
    set(gui_data.status_txt, 'String', msg);
end

%==========================================================================
% HELPERS
%==========================================================================

function im = get_slice(gui_data)
    im = volumeIdtoImage(gui_data.vol, [gui_data.curr_slice, gui_data.slice_dim]);
end

function flash_title(gui_data, msg)
    gui_data.pp.title(msg);
    drawnow;
    pause(0.35);
    gui_data.pp.title(gui_data.base_title);
end

function do_save(gui_data)
    points     = gui_data.points;      %#ok<NASGU>
    diameter_um = gui_data.diameter_um; %#ok<NASGU>
    slice_dim  = gui_data.slice_dim;   %#ok<NASGU>
    save(gui_data.save_file, 'points', 'diameter_um', 'slice_dim');
    fprintf('Saved %d points to %s\n', size(points,1), gui_data.save_file);
end

%==========================================================================
% COORDINATE CONVERSION: permuted-20um-volume indices → original voxels
%==========================================================================

function pts_orig = permuted20umToOrigVox(pts_perm, how_to_perm, ori_pxsize, registres, size_perm)
% Convert Nx3 [row col slice] indices in the permuted 20 µm registration
% volume to Nx3 original-voxel indices expected by transformPointsToAtlas.
%
% how_to_perm : signed permutation vector (sign = flip direction)
% ori_pxsize  : [row col slice] pixel size of the original sample [µm]
% registres   : registration resolution [µm], e.g. 20
% size_perm   : [dim1 dim2 dim3] size of the permuted volume

    perm_order = abs(how_to_perm);

    % Inverse permutation: inv_perm(perm_order(d)) = d
    inv_perm = zeros(1, 3);
    inv_perm(perm_order) = 1:3;

    % Step 1 – undo flips (flips were applied to permuted dims)
    perm_no_flip = pts_perm;
    for d = 1:3
        if how_to_perm(d) < 0
            perm_no_flip(:, d) = size_perm(d) + 1 - pts_perm(:, d);
        end
    end

    % Step 2 – undo permutation to recover original orientation (still 20 µm)
    N = size(pts_perm, 1);
    bvol20 = zeros(N, 3);
    for m = 1:3
        bvol20(:, m) = perm_no_flip(:, inv_perm(m));
    end

    % Step 3 – scale from 20 µm voxel indices to original voxel indices
    pts_orig = (bvol20 - 1) .* (registres ./ ori_pxsize) + 1;
end

%==========================================================================
% FIT AND VISUALISE
%==========================================================================

function fit_and_show_atlas(gui_fig)
    gui_data = guidata(gui_fig);

    if size(gui_data.points, 1) < 4
        warndlg('At least 4 edge points are needed to fit a fiber. Add more points.', ...
                'Not enough points');
        return;
    end

    fprintf('--- GRIN Lens Fitting ---\n');
    fprintf('Converting %d points to atlas space...\n', size(gui_data.points, 1));

    % Convert permuted-20µm indices → original sample voxel indices
    size_perm = [size(gui_data.vol,1), size(gui_data.vol,2), size(gui_data.vol,3)];
    pts_orig  = permuted20umToOrigVox( ...
        gui_data.points, ...
        gui_data.how_to_perm, ...
        gui_data.ori_pxsize, ...
        gui_data.registres, ...
        size_perm);

    % Transform to Allen CCF atlas space (10 µm voxels)
    atlas_pts = transformPointsToAtlas(pts_orig, 'transform_params', gui_data.trstruct);
    atlas_pts = atlas_pts(:, 1:3);

    fprintf('Atlas coords: AP [%.0f–%.0f]  DV [%.0f–%.0f]  ML [%.0f–%.0f]\n', ...
        min(atlas_pts(:,1)), max(atlas_pts(:,1)), ...
        min(atlas_pts(:,2)), max(atlas_pts(:,2)), ...
        min(atlas_pts(:,3)), max(atlas_pts(:,3)));

    % Fiber radius in 10 µm atlas voxels
    radius_vox = (gui_data.diameter_um / 2) / 10;

    fprintf('Fitting fiber circle (radius %.1f vox = %.0f µm)...\n', ...
        radius_vox, gui_data.diameter_um/2);

    [center_vox, normal_vox] = fitFiberInAtlas(atlas_pts, radius_vox);

    fprintf('Center:  [%.1f  %.1f  %.1f]  (AP DV ML in 10 µm vox)\n', center_vox);
    fprintf('Normal:  [%.3f  %.3f  %.3f]\n', normal_vox);

    % Load annotation volume
    allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
    if isempty(allen_atlas_path)
        error('annotation_10.nii.gz not found on MATLAB path.');
    end
    fprintf('Loading annotation volume...\n');
    av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));

    % Load parcellation info for region names
    parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));
    parcelinfo  = readtable(fullfile(parcel_path, 'parcellation_to_parcellation_term_membership.csv'));
    substridx   = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib]  = unique(parcelinfo.parcellation_index(substridx));
    namessub    = parcelinfo.parcellation_term_name(substridx);
    namessub    = namessub(ib);

    % Extract annotation slices at 0 : 50 : 250 µm from lens bottom
    depths_um  = 0 : 50 : 250;
    depths_vox = depths_um / 10;        % 10 µm/voxel
    Ndepths    = numel(depths_um);

    fprintf('Extracting %d annotation slices...\n', Ndepths);
    slices_av = cell(1, Ndepths);
    rvec_arr  = cell(1, Ndepths);
    for ii = 1:Ndepths
        [slices_av{ii}, rvec_arr{ii}] = extractAnnotationCircularSlice( ...
            av, center_vox, normal_vox, radius_vox, depths_vox(ii));
    end

    % Save results alongside the GUI save file
    results.center_vox  = center_vox;
    results.normal_vox  = normal_vox;
    results.radius_vox  = radius_vox;
    results.atlas_pts   = atlas_pts;
    results.depths_um   = depths_um;
    results.slices_av   = slices_av;
    results.rvec_arr    = rvec_arr;
    outfile = fullfile(gui_data.savepath, 'grin_fiber_atlas.mat');
    save(outfile, '-struct', 'results');
    fprintf('Atlas results saved to %s\n', outfile);

    % Show atlas region figure
    plotGRINAtlasImages(slices_av, rvec_arr, radius_vox, depths_um, areaidx, namessub);
end
