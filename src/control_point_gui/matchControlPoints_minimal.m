function matchControlPoints_minimal(opts)
% MATCHCONTROLPOINTS_MINIMAL Minimal control-point GUI for two volumes.
%   Shows a sample volume (left) and an atlas volume (right) side by side and
%   lets the user place matched control points. Both volumes are sliced along
%   their first dimension.
%
%   Voxel sizes are used ONLY to set the display aspect ratio of each panel,
%   never to resample the volumes.
%
%   The sample display contrast comes from per-channel intensity quantiles
%   (default [0.1 0.95]). Press Q at any time to set your own [low high]
%   quantiles; the sample panel re-normalizes without reloading the data.
%
%   Every presented half is an independent data point: the four halves of a
%   sample slice may match different atlas slices when anchored, and get
%   separate predictions when not. Halves never borrow each other's slice.
%   Per half, the atlas slice is taken from, in order:
%       1. atlas points on this half (a saved session always reopens on
%          exactly the atlas slice it was annotated on)
%       2. linear interpolation between this half's anchors (a sample slice
%          of that half carrying atlas points) that bracket this slice
%       3. the affine fit, once >=5 non-coplanar pairs exist, mapped through
%          this half's own visible region
%       4. the nearest anchor's constant offset, else the current scroll offset
%
%   Required fields of opts:
%       opts.sample      - sample volume (3D numeric), sliced along dim 1
%       opts.atlas       - atlas template volume (3D numeric), sliced along dim 1
%       opts.annotation  - atlas annotation/label volume (3D numeric)
%       opts.sample_res  - sample voxel size, scalar or [z y x]
%       opts.atlas_res   - atlas voxel size, scalar or [z y x]
%       opts.savepath    - .mat file (or folder) to save control points
%
%   Part of the LightSuite toolbox.

%==========================================================================
% 1. Options and save target
%==========================================================================
reqfields = {'sample','atlas','annotation','sample_res','atlas_res','savepath'};
for k = 1:numel(reqfields)
    if ~isfield(opts, reqfields{k})
        error('matchControlPoints_minimal:missingField', ...
            'opts is missing required field "%s".', reqfields{k});
    end
end

gui_data = struct;

[savedir, savename, saveext] = fileparts(opts.savepath);
if isempty(saveext)
    gui_data.save_path     = opts.savepath;
    gui_data.save_filename = 'control_points_minimal.mat';
else
    gui_data.save_path     = savedir;
    gui_data.save_filename = [savename saveext];
end
if ~isempty(gui_data.save_path) && ~exist(gui_data.save_path, 'dir')
    mkdir(gui_data.save_path);
end

%==========================================================================
% 2. Volumes (no resampling) and display aspect ratios
%==========================================================================
permvec         = opts.permute_sample_to_atlas;
% The sample may carry multiple channels in its 4th dimension
% (Ny x Nx x Nz x Nchannels). Permute each channel independently so the
% slicing (dim 1) axis matches the atlas, and keep the permuted raw volume
% so the display quantiles can be changed live (Q) without touching opts.
nch = size(opts.sample, 4);
rawch = cell(1, nch);
for c = 1:nch
    rawch{c} = permuteBrainVolume(single(opts.sample(:,:,:,c)), permvec);
end
gui_data.volume_raw       = cat(4, rawch{:});
gui_data.sample_quantiles = [0.1 0.95];   % [low high] display quantiles, set with Q
gui_data.volume           = requantize_sample(gui_data.volume_raw, gui_data.sample_quantiles);
gui_data.nchannels = nch;
gui_data.curr_chan = 1;   % 1-9 = single channel; 0 = RGB (first 3 channels)
gui_data.atlas  = to_uint8(opts.atlas);
gui_data.av     = resampleAnnotation(opts.annotation, opts.parcelinfo, 'structure');  
% kept for later overlay features
% Slices are taken along dim 1, so the displayed image has rows = dim 2 and
% columns = dim 3. DataAspectRatio [1 k 1] with k = res_col/res_row gives
% physically correct proportions.
sres = normalize_res(opts.sample_res(abs(permvec)));
ares = normalize_res(opts.atlas_res);
gui_data.sample_daspect = [1, sres(3)/sres(2), 1];
gui_data.atlas_daspect  = [1, ares(3)/ares(2), 1];

gui_data.nsamp  = size(gui_data.volume, 1);
gui_data.natlas = size(gui_data.atlas, 1);

%==========================================================================
% 3. Control points and slice state
%==========================================================================
% Present each sample slice 4 times, each showing ~60% of the image with one
% side blacked out: half 1=left, 2=right, 3=top, 4=bottom (labels a,b,c,d).
% chooselist rows: [sample_slice, half]; ordered slice1 a,b,c,d, slice2 ...
gui_data.nhalves = 4;
[ss, hh] = ndgrid(1:gui_data.nsamp, 1:gui_data.nhalves);
gui_data.chooselist = sortrows([ss(:), hh(:)], [1 2]);
gui_data.nentries   = size(gui_data.chooselist, 1);

% One cell entry per presented half. Points are [slice, row(dim2), col(dim3), t].
gui_data.sample_points = repmat({zeros(0,4)}, gui_data.nentries, 1);
gui_data.atlas_points  = repmat({zeros(0,4)}, gui_data.nentries, 1);

gui_data.curr_idx     = 1;                 % index into chooselist
gui_data.slice_offset = 0;                 % atlas = sample + offset (pre-anchor)
gui_data.atlas_slice  = 1;

% Affine fit state (atlas -> sample). Populated once >=5 non-coplanar pairs.
gui_data.tform    = affinetform3d(eye(4));
gui_data.has_fit  = false;
gui_data.fit_str  = '';                    % fit status line, always kept in title
gui_data.volwrap  = [];                    % annotation warped into sample space
gui_data.Rvolume  = imref3d(size(gui_data.volume, 1:3));  % sample grid
gui_data.Rmoving  = imref3d(size(gui_data.av));      % atlas grid

load_fn = fullfile(gui_data.save_path, gui_data.save_filename);
if exist(load_fn, 'file')
    old = load(load_fn);
    if isfield(old, 'sample_points') && numel(old.sample_points) == gui_data.nentries
        gui_data.sample_points = old.sample_points;
        gui_data.atlas_points  = old.atlas_points;
    end
end

%==========================================================================
% 4. GUI setup
%==========================================================================
screen_size_px = get(0, 'screensize');
gui_aspect     = 2.0;
gui_width_px   = screen_size_px(3) * 0.6;
gui_position   = [ ...
    (screen_size_px(3)-gui_width_px)/2, ...
    (screen_size_px(4)-gui_width_px/gui_aspect)/2, ...
    gui_width_px, gui_width_px/gui_aspect];

gui_fig = figure('KeyPressFcn', @keypress, ...
    'WindowScrollWheelFcn', @scroll_atlas_slice, ...
    'Toolbar', 'none', 'Menubar', 'none', 'color', 'w', ...
    'Units', 'pixels', 'Position', gui_position, ...
    'CloseRequestFcn', @close_gui);

gui_data.pp = panel();
gui_data.pp.pack('h', 2);
gui_data.pp.margin = [1 1 1 25];

controls_str = ['\bfControls: \rm\leftarrow/\rightarrow: sample slice | ' ...
    'Scroll: atlas slice | Click: add point | Backspace: delete | ' ...
    'C: clear slice | Space: overlay + annotation | 1-9/0: sample channel (0=RGB) | ' ...
    'Q: sample quantiles | S: save'];
gui_data.base_title = {'\bfControl point GUI: \rmmatch points between sample and atlas', ...
    controls_str};
gui_data.pp.title(gui_data.base_title);
gui_data.pp.fontname = 'Arial';

% --- Sample axis ---
gui_data.sample_ax = gui_data.pp(1).select();
gui_data.sample_ax.YDir = 'reverse';
gui_data.sample_ax.Colormap = gray;
hold(gui_data.sample_ax, 'on');
gui_data.sample_im_h = imagesc(get_sample_image(gui_data, 1, 1), ...
    'Parent', gui_data.sample_ax, 'ButtonDownFcn', @mouseclick_sample);
set(gui_data.sample_ax, 'DataAspectRatio', gui_data.sample_daspect);
axis(gui_data.sample_ax, 'tight');
axis(gui_data.sample_ax, 'off');

% --- Atlas axis ---
gui_data.atlas_ax = gui_data.pp(2).select();
gui_data.atlas_ax.YDir = 'reverse';
gui_data.atlas_ax.Colormap = gray;
hold(gui_data.atlas_ax, 'on');
gui_data.atlas_im_h = imagesc(volumeIdtoImage(gui_data.atlas, [1 1]), ...
    'Parent', gui_data.atlas_ax, 'ButtonDownFcn', @mouseclick_atlas);
set(gui_data.atlas_ax, 'DataAspectRatio', gui_data.atlas_daspect);
axis(gui_data.atlas_ax, 'tight');
axis(gui_data.atlas_ax, 'off');

% Warped-atlas boundary overlay (on the sample panel)
gui_data.h_overlay = scatter(gui_data.sample_ax, nan, nan,3, 'filled','MarkerFaceColor','r', ...
     'PickableParts', 'none');

% Annotation boundary overlay (on the atlas panel), toggled together with the
% sample overlay via Space.
gui_data.h_overlay_atlas = scatter(gui_data.atlas_ax, nan, nan,3, 'filled','MarkerFaceColor','c',...
     'PickableParts', 'none');

% Point marker handles
gui_data.h_pts_samp  = plot(gui_data.sample_ax, nan, nan, '.g', 'MarkerSize', 20);
gui_data.h_pts_atlas = plot(gui_data.atlas_ax,  nan, nan, '.r', 'MarkerSize', 20);
gui_data.h_text_samp  = gobjects(0);
gui_data.h_text_atlas = gobjects(0);

guidata(gui_fig, gui_data);
update_fit(gui_fig);
update_slice(gui_fig);

end

%==========================================================================
% Helpers
%==========================================================================
function res = normalize_res(res)
    res = res(:)';
    if isscalar(res), res = res([1 1 1]); end
    if numel(res) ~= 3
        error('Voxel size must be a scalar or a 3-element vector.');
    end
end

function v = to_uint8(v, quantsuse)
    if nargin < 2
        quantsuse = [0.01 0.99];
    end
    v = single(v);
    hilow = quantile(v(:), quantsuse);
    % if hi <= 0, hi = max(v(:)); end
    % if hi <= 0, hi = 1; end
    v = uint8((v - hilow(1))/range(hilow) * 255);
end

function vol = requantize_sample(raw, quants)
    % Map the permuted raw sample to a uint8 display volume, normalizing each
    % channel independently to the [low high] quantiles. Called at setup and
    % whenever the user changes the quantiles (Q).
    nch = size(raw, 4);
    vol = zeros(size(raw), 'uint8');
    for c = 1:nch
        vol(:,:,:,c) = to_uint8(raw(:,:,:,c), quants);
    end
end

function q = ask_sample_quantiles(current)
    % Prompt for the sample display quantiles. Returns [low high] on a valid
    % entry, or [] if the user cancelled or gave something unusable.
    q  = [];
    in = inputdlg({'Lower quantile (0-1):', 'Upper quantile (0-1):'}, ...
        'Sample display quantiles', [1 40], ...
        {num2str(current(1)), num2str(current(2))});
    if isempty(in), return; end
    lo = str2double(in{1});
    hi = str2double(in{2});
    if any(isnan([lo hi])) || lo < 0 || hi > 1 || lo >= hi
        warndlg('Enter two increasing values in [0, 1] with low < high.', ...
            'Invalid quantiles');
        return;
    end
    q = [lo hi];
end

function [sVals, aVals] = anchor_pairs_half(gui_data, half)
    % Anchors for a single half: sample slices (of that half) that carry atlas
    % points, with the anchor's atlas slice = median plane of those points.
    % Each half keeps its own correspondence, so halves can map to different
    % atlas planes/angles.
    sVals = [];
    aVals = [];
    idxs = find(gui_data.chooselist(:,2) == half)';
    for k = idxs
        ap = gui_data.atlas_points{k};
        if ~isempty(ap)
            sVals(end+1,1) = gui_data.chooselist(k,1);  %#ok<AGROW>
            aVals(end+1,1) = round(median(ap(:,1)));    %#ok<AGROW>
        end
    end
end

function tf = any_atlas_points(gui_data)
    tf = any(~cellfun(@isempty, gui_data.atlas_points));
end

function [d2c, d3c] = half_center(nd2, nd3, half)
    % Center of the visible ~60% region for a half, in image coords
    % (row = dim2, col = dim3).
    switch half
        case 1, d2c = nd2/2;      d3c = 0.3*nd3;   % left
        case 2, d2c = nd2/2;      d3c = 0.7*nd3;   % right
        case 3, d2c = 0.3*nd2;    d3c = nd3/2;     % top
        case 4, d2c = 0.7*nd2;    d3c = nd3/2;     % bottom
        otherwise, d2c = nd2/2;   d3c = nd3/2;
    end
end

function m = half_mask(sz, half)
    % Logical mask of the ~60% region kept for a given half; rest is blacked.
    nr = sz(1); nc = sz(2); f = 0.6;
    m = false(nr, nc);
    switch half
        case 1, m(:, 1:round(f*nc))          = true;   % left
        case 2, m(:, nc-round(f*nc)+1:nc)     = true;   % right
        case 3, m(1:round(f*nr), :)          = true;   % top
        case 4, m(nr-round(f*nr)+1:nr, :)     = true;   % bottom
    end
end

function im = get_sample_image(gui_data, s, half)
    % Displayable sample slice for the current channel selection, with the
    % complementary half blacked out. curr_chan: 1-9 = single channel (contrast
    % enhanced); 0 = truecolor RGB from the first three channels.
    raw = volumeIdtoImage(gui_data.volume, [s 1]);   % Nd2 x Nd3 (x Nchannels)
    if gui_data.curr_chan == 0
        nc3 = min(size(raw, 3), 3);
        im  = zeros([size(raw,1), size(raw,2), 3], 'uint8');
        for c = 1:nc3
            % im(:,:,c) = adapthisteq(raw(:,:,c));
            im(:,:,c) = raw(:,:,c);
        end
    else
        c  = min(gui_data.curr_chan, size(raw,3));
        % im = adapthisteq(raw(:,:,c));
        im = raw(:,:,c);
    end
    im = apply_half_mask(im, half);
end

function im = apply_half_mask(im, half)
    % Blank the complementary region for every channel of im.
    m = half_mask([size(im,1) size(im,2)], half);
    for c = 1:size(im,3)
        ch = im(:,:,c);
        ch(~m) = 0;
        im(:,:,c) = ch;
    end
end

function draw_atlas_annotation(gui_data)
    % Region boundaries of the annotation volume at the current atlas slice,
    % over the full slice. Visibility is controlled via Space.
    a    = min(max(gui_data.atlas_slice, 1), gui_data.natlas);
    ann  = volumeIdtoImage(gui_data.av, [a 1]);
    % ann  = medfilt2(ann, [5 5]);
    bnd  = round(conv2(single(ann), ones(3)./9, 'same')) ~= ann;
    % bnd  = imgradient(ann)~=0;
    [row, col] = ind2sub(size(ann), find(bnd));
    set(gui_data.h_overlay_atlas, 'XData', col, 'YData', row);
end

function a = atlas_slice_for(gui_data, idx)
    % Atlas slice for a presented half (chooselist entry). Each half is its own
    % data point: its own points are authoritative, so an annotated half always
    % comes back on the atlas slice it was annotated on, and an unanchored half
    % is predicted from its own correspondence alone. Halves never borrow each
    % other's slice, so they can sit on different atlas planes/angles.
    s    = gui_data.chooselist(idx, 1);
    half = gui_data.chooselist(idx, 2);

    % 1. This half's own points win (explicit user correspondence).
    ap = gui_data.atlas_points{idx};
    if ~isempty(ap)
        a = clamp_slice(median(ap(:,1)), gui_data);
        return;
    end

    % 2. Otherwise interpolate THIS half's own anchors, where they bracket the
    %    slice. Only interpolation is used: extrapolating a line off a noisy
    %    tail of anchors can march backwards through the atlas.
    [sVals, aVals] = anchor_pairs_half(gui_data, half);
    a = interp_anchors(sVals, aVals, s);
    if ~isnan(a)
        a = clamp_slice(a, gui_data);
        return;
    end

    % 3. With a fit, map THIS half's visible-region center through the inverse
    %    transform, so a tilted plane yields a different slice per half.
    if gui_data.has_fit
        [d2c, d3c] = half_center(size(gui_data.volume,2), size(gui_data.volume,3), half);
        c = [d2c, s, d3c];   % intrinsic [x y z] = [dim2, dim1, dim3]
        p = gui_data.tform.transformPointsInverse(c);
        a = clamp_slice(p(2), gui_data);   % atlas dim1 (y-intrinsic)
        return;
    end

    % 4. Off the end of this half's anchors with no fit to lean on: hold the
    %    nearest anchor's offset (slope 1). With no anchors at all the atlas
    %    just follows the current scroll offset.
    if ~isempty(sVals)
        [~, inear] = min(abs(sVals - s));
        a = clamp_slice(aVals(inear) + (s - sVals(inear)), gui_data);
    else
        a = clamp_slice(s + gui_data.slice_offset, gui_data);
    end
end

function a = clamp_slice(a, gui_data)
    a = min(max(round(a), 1), gui_data.natlas);
end

function a = interp_anchors(sVals, aVals, s)
    % Anchor interpolation, NaN unless the anchors actually bracket s.
    if numel(sVals) >= 2 && s >= min(sVals) && s <= max(sVals)
        a = interp1(sVals, aVals, s, 'linear');
    else
        a = nan;
    end
end

%==========================================================================
% Navigation / drawing
%==========================================================================
function update_slice(gui_fig)
    gui_data = guidata(gui_fig);
    idx  = gui_data.curr_idx;
    s    = gui_data.chooselist(idx, 1);
    half = gui_data.chooselist(idx, 2);

    % Sample image (selected channel or RGB), complementary region blacked out.
    set(gui_data.sample_im_h, 'CData', get_sample_image(gui_data, s, half));
    labels = 'abcd';
    title(gui_data.sample_ax, sprintf('Sample slice %d/%d (%c)', ...
        s, gui_data.nsamp, labels(half)), 'FontSize', 12);

    % Atlas slice follows this half's own anchor / fit (see atlas_slice_for).
    gui_data.atlas_slice = atlas_slice_for(gui_data, idx);

    guidata(gui_fig, gui_data);
    draw_markers(gui_data, 'sample');
    draw_overlay(gui_fig);
    update_atlas_slice(gui_fig);
end

function update_atlas_slice(gui_fig)
    gui_data = guidata(gui_fig);
    a = min(max(gui_data.atlas_slice, 1), gui_data.natlas);
    gui_data.atlas_slice = a;

    % The atlas is always shown in full (no half blackout, unlike the sample).
    atlas_im = adapthisteq(volumeIdtoImage(gui_data.atlas, [a 1]));
    % atlas_im = volumeIdtoImage(gui_data.atlas, [a 1]);
    set(gui_data.atlas_im_h, 'CData', atlas_im);
    title(gui_data.atlas_ax, sprintf('Atlas slice %d/%d', a, gui_data.natlas), ...
        'FontSize', 12);

    guidata(gui_fig, gui_data);
    draw_atlas_annotation(gui_data);
    draw_markers(gui_data, 'atlas');
end

function draw_markers(gui_data, which)
    % Rows = dim 2 (y), columns = dim 3 (x).
    if strcmp(which, 'sample')
        h_plot   = gui_data.h_pts_samp;
        ax       = gui_data.sample_ax;
        pts      = gui_data.sample_points{gui_data.curr_idx};
        old_text = gui_data.h_text_samp;
    else
        h_plot   = gui_data.h_pts_atlas;
        ax       = gui_data.atlas_ax;
        % Only show atlas points that live on the currently shown atlas slice.
        allpts   = gui_data.atlas_points{gui_data.curr_idx};
        if ~isempty(allpts)
            pts = allpts(round(allpts(:,1)) == gui_data.atlas_slice, :);
        else
            pts = allpts;
        end
        old_text = gui_data.h_text_atlas;
    end

    if ~isempty(pts)
        set(h_plot, 'XData', pts(:,3), 'YData', pts(:,2));
    else
        set(h_plot, 'XData', nan, 'YData', nan);
    end

    delete(old_text(isvalid(old_text)));
    new_text = gobjects(size(pts,1), 1);
    for i = 1:size(pts,1)
        new_text(i) = text(ax, pts(i,3), pts(i,2), num2str(i), ...
            'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
            'PickableParts', 'none');
    end

    if strcmp(which, 'sample')
        gui_data.h_text_samp = new_text;
    else
        gui_data.h_text_atlas = new_text;
    end
    guidata(gui_data.pp.figure, gui_data);
end

function scroll_atlas_slice(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    gui_data.atlas_slice = min(max( ...
        gui_data.atlas_slice + eventdata.VerticalScrollCount, 1), gui_data.natlas);
    % With no committed anchors anywhere, remember the scroll as a running offset.
    if ~any_atlas_points(gui_data)
        s = gui_data.chooselist(gui_data.curr_idx, 1);
        gui_data.slice_offset = gui_data.atlas_slice - s;
    end
    guidata(gui_fig, gui_data);
    update_atlas_slice(gui_fig);
end

%==========================================================================
% Clicks
%==========================================================================
function mouseclick_sample(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    idx = gui_data.curr_idx;
    s   = gui_data.chooselist(idx, 1);
    p = eventdata.IntersectionPoint;   % [x(col) y(row)]
    cpt = [s, p(2), p(1), now];
    gui_data.sample_points{idx} = [gui_data.sample_points{idx}; cpt];
    guidata(gui_fig, gui_data);
    draw_markers(gui_data, 'sample');
    update_fit(gui_fig);
end

function mouseclick_atlas(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    idx = gui_data.curr_idx;
    p = eventdata.IntersectionPoint;   % [x(col) y(row)]
    cpt = [gui_data.atlas_slice, p(2), p(1), now];
    gui_data.atlas_points{idx} = [gui_data.atlas_points{idx}; cpt];
    guidata(gui_fig, gui_data);
    draw_markers(gui_data, 'atlas');
    update_fit(gui_fig);
end

%==========================================================================
% Keyboard
%==========================================================================
function keypress(gui_fig, eventdata)
    gui_data = guidata(gui_fig);
    needs_update = false;
    needs_refit  = false;

    switch eventdata.Key
        case 'leftarrow'
            gui_data.curr_idx = max(gui_data.curr_idx - 1, 1);
            needs_update = true;
        case 'rightarrow'
            gui_data.curr_idx = min(gui_data.curr_idx + 1, gui_data.nentries);
            needs_update = true;
        case 'space'
            vis = get(gui_data.h_overlay, 'Visible');
            newvis = cell2mat(setdiff({'on','off'}, vis));
            set(gui_data.h_overlay, 'Visible', newvis);
            set(gui_data.h_overlay_atlas, 'Visible', newvis);
        case {'1','2','3','4','5','6','7','8','9'}
            ch = str2double(eventdata.Key);
            if ch <= gui_data.nchannels
                gui_data.curr_chan = ch;
                needs_update = true;
            end
        case '0'
            if gui_data.nchannels >= 3
                gui_data.curr_chan = 0;   % RGB from first three channels
                needs_update = true;
            end
        case 'return'
            in = inputdlg(sprintf('Go to sample slice (1-%d):', gui_data.nsamp));
            if ~isempty(in)
                v = round(str2double(in{1}));
                if ~isnan(v) && v >= 1 && v <= gui_data.nsamp
                    gui_data.curr_idx = (v-1)*gui_data.nhalves + 1;  % first half
                    needs_update = true;
                end
            end
        case 'q'
            q = ask_sample_quantiles(gui_data.sample_quantiles);
            if ~isempty(q)
                gui_data.sample_quantiles = q;
                gui_data.volume = requantize_sample(gui_data.volume_raw, q);
                needs_update = true;
            end
        case 'c'
            gui_data.sample_points{gui_data.curr_idx} = zeros(0,4);
            gui_data.atlas_points{gui_data.curr_idx}  = zeros(0,4);
            needs_update = true;
            needs_refit  = true;
        case 'backspace'
            sp = gui_data.sample_points{gui_data.curr_idx};
            ap = gui_data.atlas_points{gui_data.curr_idx};
            t_s = -inf; t_a = -inf;
            if ~isempty(sp), t_s = sp(end,4); end
            if ~isempty(ap), t_a = ap(end,4); end
            if t_s > t_a
                gui_data.sample_points{gui_data.curr_idx}(end,:) = [];
            elseif t_a > -inf
                gui_data.atlas_points{gui_data.curr_idx}(end,:) = [];
            end
            needs_update = true;
            needs_refit  = true;
        case 's'
            save_fn = savedata(gui_data);
            gui_data.pp.title({'\bfSAVED!'});
            pause(0.2);
            set_title(gui_data);   % restore controls + fit status line
            disp(['Saved ' save_fn]);
    end

    if needs_update
        guidata(gui_fig, gui_data);
        if needs_refit
            update_fit(gui_fig);
        end
        update_slice(gui_fig);
    end
end

%==========================================================================
% Affine fit (atlas -> sample)
%==========================================================================
function [src, tgt, n] = matched_pairs(gui_data)
    % Pair sample/atlas points by click order within each slice. Points are
    % stored [d1 d2 d3 t]; reorder to intrinsic [x y z] = [d2 d1 d3] for the
    % transform / imwarp conventions.
    src = zeros(0,3);   % atlas
    tgt = zeros(0,3);   % sample
    for i = 1:gui_data.nentries
        sp = gui_data.sample_points{i};
        ap = gui_data.atlas_points{i};
        k = min(size(sp,1), size(ap,1));
        if k > 0
            tgt = [tgt; sp(1:k, [2 1 3])]; %#ok<AGROW>
            src = [src; ap(1:k, [2 1 3])]; %#ok<AGROW>
        end
    end
    n = size(src, 1);
end

function tf = is_noncoplanar(pts)
    % True when the points span all three dimensions (needed for affine).
    if size(pts,1) < 4, tf = false; return; end
    s = svd(pts - mean(pts,1));
    tf = s(3) > 1e-6 * s(1);
end

function update_fit(gui_fig)
    gui_data = guidata(gui_fig);
    [src, tgt, n] = matched_pairs(gui_data);

    if n >= 5 && is_noncoplanar(src)
        if n >=16
            [tform, mse] = fitAffineTrans3D(src, tgt);
        else
            [tform, mse] = fitSimilarityTrans3D(src, tgt);
        end
        gui_data.tform   = tform;
        gui_data.has_fit = true;
        gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, ...
            'nearest', 'OutputView', gui_data.Rvolume);
        gui_data.fit_str = sprintf('\\bfAffine fit: \\rmMSE %.2f | N %d', mse, n);
    else
        gui_data.has_fit = false;
        gui_data.volwrap = [];
        gui_data.fit_str = sprintf( ...
            '\\bfAffine fit: \\rmneed >= 5 non-coplanar pairs (have %d)', n);
    end
    set_title(gui_data);

    guidata(gui_fig, gui_data);
    draw_overlay(gui_fig);
end

function set_title(gui_data)
    % Title = base controls + persistent fit status (points and MSE), so the
    % fit line stays visible through slice changes, saves, etc.
    gui_data.pp.title([gui_data.base_title, {gui_data.fit_str}]);
end

function draw_overlay(gui_fig)
    gui_data = guidata(gui_fig);
    if isempty(gui_data.volwrap)
        set(gui_data.h_overlay, 'XData', nan, 'YData', nan);
        return;
    end
    idx  = gui_data.curr_idx;
    s    = gui_data.chooselist(idx, 1);
    half = gui_data.chooselist(idx, 2);
    slicewarp = volumeIdtoImage(gui_data.volwrap, [s 1]);
    bnd = round(conv2(single(slicewarp), ones(3)./9, 'same')) ~= slicewarp;
    bnd = bnd & half_mask(size(slicewarp), half);   % keep only the shown region
    [row, col] = ind2sub(size(slicewarp), find(bnd));
    set(gui_data.h_overlay, 'XData', col, 'YData', row);
end

%==========================================================================
% Save / close
%==========================================================================
function save_fn = savedata(gui_data)
    sample_points = gui_data.sample_points;   %#ok<NASGU>
    atlas_points  = gui_data.atlas_points;    %#ok<NASGU>
    % Affine mapping atlas -> sample (intrinsic [x y z]); identity if unfit.
    atlas2sample_tform = gui_data.tform.A;    %#ok<NASGU>
    has_fit = gui_data.has_fit;               %#ok<NASGU>
    save_fn = fullfile(gui_data.save_path, gui_data.save_filename);
    save(save_fn, 'sample_points', 'atlas_points', 'atlas2sample_tform', 'has_fit');
end

function close_gui(gui_fig, ~)
    gui_data = guidata(gui_fig);
    switch questdlg('Save changes?', 'Confirm exit', 'Yes', 'No', 'Cancel', 'Yes')
        case 'Yes'
            save_fn = savedata(gui_data);
            disp(['Saved ' save_fn]);
            delete(gui_fig);
        case 'No'
            delete(gui_fig);
    end
end
