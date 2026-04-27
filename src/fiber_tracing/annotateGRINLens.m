function annotateGRINLens(savepath, varargin)
% ANNOTATEGRINLENS  GUI for annotating optical fiber / GRIN lens positions
%   in a LightSuite registered volume.
%
%   Displays the 20 µm registration volume in coronal view.  The user
%   navigates slices and clicks the fiber edges.  Up to 9 fibers can be
%   annotated simultaneously, each in a separate numbered group.  Pressing F
%   transforms all groups to Allen-CCF atlas space, fits a circle per fiber,
%   and opens one atlas-region figure per fiber (depths 0–250 µm).
%
%   Usage:
%       annotateGRINLens(savepath)
%       annotateGRINLens(savepath, 'Diameter', 500)
%
%   Required:
%       savepath   – LightSuite folder with regopts.mat + transform_params.mat
%
%   Optional:
%       Diameter   – fiber diameter in µm  (default 500)
%
%   Keyboard controls:
%       ←  / →  or scroll   navigate coronal slices
%       1 – 9               select active fiber group (colour-coded)
%       Click               add edge point to active fiber
%       Backspace           delete last point from active fiber
%       C                   clear active fiber's points on current slice
%       S                   save all points to grin_fiber_points.mat
%       F                   fit every non-empty fiber + show atlas images
%       Enter               jump to a specific slice number

p = inputParser;
addRequired(p,  'savepath');
addParameter(p, 'Diameter', 500, @isnumeric);
parse(p, savepath, varargin{:});

% 9 distinct group colours (used for point overlays and status labels)
GROUP_COLORS = [ ...
    1.00 0.20 0.20;   % 1  red
    0.15 0.82 0.15;   % 2  green
    0.28 0.50 1.00;   % 3  blue
    1.00 0.55 0.00;   % 4  orange
    0.82 0.15 0.82;   % 5  magenta
    0.00 0.78 0.88;   % 6  cyan
    0.80 0.75 0.00;   % 7  yellow
    0.55 0.00 0.85;   % 8  purple
    0.00 0.65 0.45];  % 9  teal

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

% Locate registration volume – resolve path relative to savepath
[~, currname, currext] = fileparts(regopts.regvolpath);
regvolpath = fullfile(savepath, [currname currext]);
if ~exist(regvolpath, 'file')
    error('annotateGRINLens: registration volume not found at %s', regvolpath);
end

fprintf('Loading registration volume from %s ...\n', regvolpath);
volraw = readDownStack(regvolpath);

% Permute to atlas orientation (same convention as the registration pipeline)
how_to_perm = trstruct.how_to_perm;
volperm     = permuteBrainVolume(volraw, how_to_perm);

% Normalise to uint8 for display
irand   = randperm(numel(volperm), min(numel(volperm), 1e5));
factv   = 255 / single(quantile(single(volperm(irand)), 0.999));
voldisp = uint8(min(single(volperm) * factv, 255));

fprintf('Volume size (permuted): %d x %d x %d\n', size(voldisp,1), size(voldisp,2), size(voldisp,3));

%----------------------------------------------------------------------
% 2. Initialise GUI state
%----------------------------------------------------------------------
NGROUPS   = 9;
SLICE_DIM = 1;           % coronal: navigate along the AP dimension

gui_data = struct;
gui_data.savepath      = savepath;
gui_data.trstruct      = trstruct;
gui_data.regopts       = regopts;
gui_data.vol           = voldisp;
gui_data.how_to_perm   = how_to_perm;
gui_data.diameter_um   = p.Results.Diameter;
gui_data.registres     = regopts.registres;
gui_data.ori_pxsize    = trstruct.ori_pxsize;
gui_data.slice_dim     = SLICE_DIM;
gui_data.curr_slice    = round(size(voldisp, SLICE_DIM) / 2);
gui_data.active_group  = 1;
gui_data.all_points    = repmat({zeros(0,3)}, 1, NGROUPS);  % {1×9} of Nx3
gui_data.group_colors  = GROUP_COLORS;
gui_data.save_file     = fullfile(savepath, 'grin_fiber_points.mat');

% Re-load previously saved points if present
if exist(gui_data.save_file, 'file')
    sv = load(gui_data.save_file, 'all_points', 'active_group');
    if isfield(sv, 'all_points') && iscell(sv.all_points)
        % Pad to NGROUPS in case the file was saved with fewer groups
        for g = 1:min(numel(sv.all_points), NGROUPS)
            gui_data.all_points{g} = sv.all_points{g};
        end
        total = sum(cellfun(@(x) size(x,1), gui_data.all_points));
        fprintf('Loaded %d points across all groups.\n', total);
    end
    if isfield(sv, 'active_group');   gui_data.active_group = sv.active_group; end
    if ~isempty(sv.all_points{1});
        gui_data.curr_slice = round(median(sv.all_points{1}(:,1)));
    end
    
end

%----------------------------------------------------------------------
% 3. Build figure and layout
%----------------------------------------------------------------------
screen_sz = get(0, 'screensize');
fig_w = screen_sz(3) * 0.55;
fig_h = fig_w * 0.72;
gui_fig = figure( ...
    'Name',                 'GRIN Lens Annotation', ...
    'KeyPressFcn',          @keypress, ...
    'WindowScrollWheelFcn', @scroll_slice, ...
    'Toolbar',              'none', ...
    'Menubar',              'none', ...
    'Color',                'w', ...
    'Units',                'pixels', ...
    'Position',             [(screen_sz(3)-fig_w)/2, (screen_sz(4)-fig_h)/2, fig_w, fig_h], ...
    'CloseRequestFcn',      @close_gui);

gui_data.pp = panel();
gui_data.pp.pack('v', {0.91, 0.09});
gui_data.pp.margin = [2 2 2 30];

gui_data.base_title = { ...
    ['\bfGRIN Lens Annotation  (diameter: ' num2str(gui_data.diameter_um) ' µm)'], ...
    ['\bfControls:\rm  ←/→ or scroll: slice  |  Click: add point  |  ' ...
     'Backspace: delete  |  C: clear slice  |  1–9: select fiber  |  ' ...
     'F: fit atlas  |  S: save  |  Enter: jump']};
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

% One coloured plot handle per group (all initialised to NaN)
gui_data.pts_h = gobjects(1, NGROUPS);
for g = 1:NGROUPS
    gui_data.pts_h(g) = plot(ax, nan, nan, '.', ...
        'Color',          GROUP_COLORS(g,:), ...
        'MarkerSize',     16, ...
        'PickableParts',  'none', ...
        'DisplayName',    sprintf('Fiber %d', g));
end

% Status line
stat_ax = gui_data.pp(2).select();
axis(stat_ax, 'off');
gui_data.status_txt = text(stat_ax, 0.01, 0.5, '', ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'middle', ...
    'Interpreter', 'none');

guidata(gui_fig, gui_data);
update_display(gui_fig);

end  % annotateGRINLens

%==========================================================================
% CALLBACKS
%==========================================================================

function keypress(gui_fig, event)
    gui_data     = guidata(gui_fig);
    needs_update = true;

    % Number keys 1-9: select active fiber group
    num_keys = {'1','2','3','4','5','6','7','8','9', ...
                'numpad1','numpad2','numpad3','numpad4','numpad5', ...
                'numpad6','numpad7','numpad8','numpad9'};
    key_vals = [1 2 3 4 5 6 7 8 9  1 2 3 4 5 6 7 8 9];
    grp_idx  = find(strcmp(event.Key, num_keys), 1);
    if ~isempty(grp_idx)
        gui_data.active_group = key_vals(grp_idx);
        guidata(gui_fig, gui_data);
        update_display(gui_fig);
        return;
    end

    switch event.Key
        case 'leftarrow'
            gui_data.curr_slice = max(1, gui_data.curr_slice - 1);

        case 'rightarrow'
            max_sl = size(gui_data.vol, gui_data.slice_dim);
            gui_data.curr_slice = min(max_sl, gui_data.curr_slice + 1);

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
            g = gui_data.active_group;
            if ~isempty(gui_data.all_points{g})
                gui_data.all_points{g}(end, :) = [];
            end

        case 'c'
            g    = gui_data.active_group;
            idim = gui_data.slice_dim;
            keep = gui_data.all_points{g}(:, idim) ~= gui_data.curr_slice;
            gui_data.all_points{g} = gui_data.all_points{g}(keep, :);

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
    max_sl   = size(gui_data.vol, gui_data.slice_dim);
    gui_data.curr_slice = max(1, min(max_sl, ...
        gui_data.curr_slice + event.VerticalScrollCount));
    guidata(gui_fig, gui_data);
    update_display(gui_fig);
end

function mouseclick(src, event)
    if ~strcmp(get(ancestor(src,'figure'), 'SelectionType'), 'normal')
        return;
    end
    gui_data  = guidata(src);
    idim      = gui_data.slice_dim;  % = 1 (coronal)
    alldims   = 1:3;
    free_dims = alldims(alldims ~= idim);   % [2, 3]
    pt3d      = zeros(1, 3);

    % IntersectionPoint is [col, row, depth] in axes data coords
    xy = event.IntersectionPoint(1:2);
    pt3d(free_dims(1)) = xy(2);    % row → dim 2
    pt3d(free_dims(2)) = xy(1);    % col → dim 3
    pt3d(idim)         = gui_data.curr_slice;

    g = gui_data.active_group;
    gui_data.all_points{g} = [gui_data.all_points{g}; pt3d];
    guidata(src, gui_data);
    update_display(ancestor(src, 'figure'));
end

function close_gui(gui_fig, ~)
    gui_data = guidata(gui_fig);
    total_pts = sum(cellfun(@(x) size(x,1), gui_data.all_points));
    if total_pts > 0
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
    gui_data  = guidata(gui_fig);
    idim      = gui_data.slice_dim;
    alldims   = 1:3;
    free_dims = alldims(alldims ~= idim);   % [2, 3]

    % Slice image
    curr_im = get_slice(gui_data);
    try; curr_im = adapthisteq(curr_im); catch; end
    set(gui_data.im_h, 'CData', curr_im);
    set(gui_data.ax, 'XLim', [0.5, size(curr_im,2)+0.5], ...
                     'YLim', [0.5, size(curr_im,1)+0.5]);

    % Points per group on current slice
    for g = 1:numel(gui_data.all_points)
        pts_g    = gui_data.all_points{g};
        on_slice = pts_g(:, idim) == gui_data.curr_slice;
        pts_vis  = pts_g(on_slice, :);
        if ~isempty(pts_vis)
            set(gui_data.pts_h(g), ...
                'XData', pts_vis(:, free_dims(2)), ...
                'YData', pts_vis(:, free_dims(1)));
        else
            set(gui_data.pts_h(g), 'XData', nan, 'YData', nan);
        end
    end

    % Highlight active group with a larger marker
    NGROUPS = numel(gui_data.all_points);
    for g = 1:NGROUPS
        sz = 14;
        if g == gui_data.active_group; sz = 20; end
        set(gui_data.pts_h(g), 'MarkerSize', sz);
    end

    % Status line
    max_sl    = size(gui_data.vol, idim);
    grp_cnts  = cellfun(@(x) size(x,1), gui_data.all_points);
    cnt_str   = strjoin(arrayfun(@(g) sprintf('%d:%d', g, grp_cnts(g)), ...
                        1:NGROUPS, 'UniformOutput', false), '  ');
    gc        = gui_data.group_colors(gui_data.active_group, :);
    col_name  = group_color_name(gui_data.active_group);
    msg = sprintf('Coronal  %d / %d  |  Active: Fiber %d (%s)  |  Pts  %s', ...
        gui_data.curr_slice, max_sl, gui_data.active_group, col_name, cnt_str);
    set(gui_data.status_txt, 'String', msg, ...
        'Color', gc);
end

%==========================================================================
% HELPERS
%==========================================================================

function im = get_slice(gui_data)
    im = volumeIdtoImage(gui_data.vol, [gui_data.curr_slice, gui_data.slice_dim]);
end

function name = group_color_name(g)
    names = {'red','green','blue','orange','magenta','cyan','yellow','purple','teal'};
    name  = names{g};
end

function flash_title(gui_data, msg)
    gui_data.pp.title(msg);
    drawnow;
    pause(0.35);
    gui_data.pp.title(gui_data.base_title);
end

function do_save(gui_data)
    all_points  = gui_data.all_points;   %#ok<NASGU>
    diameter_um = gui_data.diameter_um;  %#ok<NASGU>
    active_group = gui_data.active_group; %#ok<NASGU>
    save(gui_data.save_file, 'all_points', 'diameter_um', 'active_group');
    total = sum(cellfun(@(x) size(x,1), gui_data.all_points));
    fprintf('Saved %d total points to %s\n', total, gui_data.save_file);
end

%==========================================================================
% FIT AND VISUALISE (all non-empty fiber groups)
%==========================================================================

function fit_and_show_atlas(gui_fig)
    gui_data = guidata(gui_fig);

    % Check that at least one group has enough points
    grp_cnts   = cellfun(@(x) size(x,1), gui_data.all_points);
    valid_grps = find(grp_cnts >= 4);
    if isempty(valid_grps)
        warndlg('Each fiber needs at least 4 edge points. Add more points and try again.', ...
                'Not enough points');
        return;
    end

    fprintf('--- GRIN Lens Fitting (%d fiber(s)) ---\n', numel(valid_grps));

    % Load atlas resources once (shared across all fibers)
    allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
    if isempty(allen_atlas_path)
        error('annotation_10.nii.gz not found on MATLAB path.');
    end
    fprintf('Loading annotation volume...\n');
    av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));

    parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));
    parcelinfo  = readtable(fullfile(parcel_path, 'parcellation_to_parcellation_term_membership.csv'));
    substridx   = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib] = unique(parcelinfo.parcellation_index(substridx));
    namessub    = parcelinfo.parcellation_term_acronym(substridx);
    namessub    = namessub(ib);

    colorsub      = [parcelinfo.red(substridx), parcelinfo.green(substridx), parcelinfo.blue(substridx)];
    colorsub     = colorsub(ib, :);

    depths_um  = 0 : 100 : 400;
    depths_vox = depths_um / 10;
    Ndepths    = numel(depths_um);

    radius_vox = (gui_data.diameter_um / 2) / 10;
    all_results = cell(1, numel(valid_grps));
 
    regsize_mm = gui_data.registres * 1e-3;
    [~, trname,trext] = fileparts(gui_data.trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    bsplinepath = fullfile(gui_data.savepath, [trname trext]);
    
    Dfield = transformix([], bsplinepath);
    Dfield = permute(Dfield, [2 3 4 1]) / regsize_mm; 
    [Sx, Sy, Sz, ~] = size(Dfield);
    
    Xgv = 1:Sx; Ygv = 1:Sy; Zgv = 1:Sz;

    for ki = 1:numel(valid_grps)
        g    = valid_grps(ki);
        ptscurr  = gui_data.all_points{g};
        fprintf('\n[Fiber %d]  %d points\n', g, size(ptscurr,1));

    
        %%

        absperm = abs(gui_data.trstruct.how_to_perm);
        facmult = gui_data.registres./gui_data.trstruct.ori_pxsize;
        facmult = facmult(absperm);
        ptsuse = ptscurr .* facmult.* gui_data.trstruct.ori_pxsize(absperm) * 1e-3; 
        pts    = ptsuse(:, [2 1 3]) / regsize_mm;

      
        dx = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,1), pts(:,1), pts(:,2), pts(:,3), 'linear');
        dy = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,2), pts(:,1), pts(:,2), pts(:,3), 'linear');
        dz = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,3), pts(:,1), pts(:,2), pts(:,3), 'linear');
        
        interpolated_displacements = -[dx, dy, dz];
        interpolated_displacements(isnan(interpolated_displacements)) = 0;
        
        atlas_pts = gui_data.trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward(pts + interpolated_displacements);

        atlas_pts = atlas_pts(:, 1:3);
        % %%
        % figure;
        % plotBrainGrid; hold on;
        % scatter3(atlas_pts(:,2),atlas_pts(:,3),atlas_pts(:,1))
        %%

     
        fprintf('  Atlas: AP [%.0f–%.0f]  DV [%.0f–%.0f]  ML [%.0f–%.0f]\n', ...
            min(atlas_pts(:,1)), max(atlas_pts(:,1)), ...
            min(atlas_pts(:,2)), max(atlas_pts(:,2)), ...
            min(atlas_pts(:,3)), max(atlas_pts(:,3)));

        [center_vox, normal_vox] = fitFiberInAtlas(atlas_pts(:, [2 1 3]), radius_vox);
        fprintf('  Centre: [%.1f  %.1f  %.1f]  Normal: [%.3f  %.3f  %.3f]\n', ...
            center_vox, normal_vox);

        % Extract annotation slices
        slices_av = cell(1, Ndepths);
        rvec_arr  = cell(1, Ndepths);
        for ii = 1:Ndepths
            [slices_av{ii}, rvec_arr{ii}] = extractAnnotationCircularSlice( ...
                av, center_vox, normal_vox, radius_vox, depths_vox(ii));
        end

        % Save per-fiber results
        res.center_vox = center_vox;
        res.normal_vox = normal_vox;
        res.radius_vox = radius_vox;
        res.atlas_pts  = atlas_pts;
        res.depths_um  = depths_um;
        res.slices_av  = slices_av;
        res.rvec_arr   = rvec_arr;
        res.atlas_pts  = atlas_pts;
        all_results{ki} = res;

        outfile = fullfile(gui_data.savepath, sprintf('grin_fiber%d_atlas.mat', g));
        save(outfile, '-struct', 'res');
        fprintf('  Results saved to %s\n', outfile);
    end

        
        % Show atlas figure for this fiber
        plotGRINAtlasImages(all_results, areaidx, namessub, colorsub);
end