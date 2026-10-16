function [intensity_results, channames] = extractGRINRawIntensitySlices(savepath, all_results, varargin)
% EXTRACTGRINRAWINTENSITYSLICES  Sample the GRIN cross-sections from the RAW,
%   full-resolution light-sheet data instead of the registered 20 µm volume.
%
%   [INTENSITY_RESULTS, CHANNAMES] = extractGRINRawIntensitySlices(SAVEPATH, ALL_RESULTS)
%
%   Same outputs as extractGRINIntensitySlices, but the cross-section grid is
%   at the native pixel size, so individual cells stay visible. The planes are
%   the same as the atlas slices (same centre, axis, basis and extent); only
%   the grid is finer.
%
%   How it stays cheap:
%     1. Every grid point (all depths of a fiber) is mapped atlas -> raw pixel
%        with atlasPointsToSample (exact inverse of the cell transform,
%        b-spline included). The displacement field is computed once.
%     2. Points are bucketed by raw plane. The raw stack is streamed once per
%        fiber: each plane that some point needs is read once, cropped to the
%        rows/cols those points span (a PixelRegion read for planeperfile
%        data), bilinearly sampled, and accumulated with its linear z weight.
%        Nothing beyond one cropped plane is held in memory.
%
%   Name-value options
%   ------------------
%   'PixelSize'  ([])     grid spacing in µm. Default: the finest raw pixel
%                         size, min(ori_pxsize).
%   'Channels'   ([])     raw channels to sample. Default: all (opts.Nchans).
%   'DataFolder' ('')     raw data folder, if it moved since processing.
%   'UseBspline' (true)   false = affine-only inverse (faster, less exact).
%   'Annotation' ([])     annotation_10 volume; when given, the parcellation
%                         is also sampled on the same fine grid (.slices_av).
%
%   Outputs
%   -------
%   INTENSITY_RESULTS – 1×Nfibers cell; each struct has slices_int (1×Ndepths
%       cell of [Ngrid×Ngrid×Nchan] single), median_intensity, mean_intensity
%       ([Ndepths×Nchan]), rvec_arr (grid vectors in 10 µm atlas voxels, so
%       they overlay the atlas slices directly), depths_um, radius_vox,
%       pxsize_um, and slices_av (annotation on the same grid, {} if no
%       'Annotation' was passed). The grid spans ±3 radii (grinSliceGrid).
%   CHANNAMES – 1×Nchan channel names.

    p = inputParser;
    addParameter(p, 'PixelSize',  [],    @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'Channels',   [],    @isnumeric);
    addParameter(p, 'DataFolder', '',    @(x) ischar(x) || isstring(x));
    addParameter(p, 'UseBspline', true,  @(x) islogical(x) || isscalar(x));
    addParameter(p, 'Annotation', [],    @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    %----------------------------------------------------------------------
    % 1. Registration options, transform and raw data
    %----------------------------------------------------------------------
    ro       = load(fullfile(savepath, 'regopts.mat'), 'opts');
    opts     = ro.opts;
    trstruct = load(fullfile(savepath, 'transform_params.mat'));

    % resolve the b-spline file relative to savepath (robust to moved folders)
    [~, bn, be] = fileparts(trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    bpath = fullfile(savepath, [bn be]);
    if exist(bpath, 'file')
        trstruct.tform_bspline_samp20um_to_atlas_20um_px = bpath;
    end

    tfiles = resolveRawDataFiles(opts, opt.DataFolder);

    Nchans = getOr(opts, 'Nchans', 1);
    chans  = opt.Channels;
    if isempty(chans); chans = 1:Nchans; end
    channames = cell(1, numel(chans));
    for k = 1:numel(chans)
        channames{k} = sprintf('ch%d', chans(k));
        if isfield(opts, 'channames') && numel(opts.channames) >= chans(k) ...
                && ~isempty(opts.channames{chans(k)})
            channames{k} = char(opts.channames{chans(k)});
        end
    end

    pxsize_um = opt.PixelSize;
    if isempty(pxsize_um); pxsize_um = min(trstruct.ori_pxsize); end
    step_vox = pxsize_um / 10;            % grid step in 10 µm atlas voxels

    Ny = trstruct.ori_size(1);  Nx = trstruct.ori_size(2);  Nz = trstruct.ori_size(3);

    fprintf('Sampling raw data at %.2f µm (raw px [%s] µm, %d channel(s)):\n', ...
        pxsize_um, num2str(trstruct.ori_pxsize, '%.2f '), numel(chans));

    %----------------------------------------------------------------------
    % 2. Per fiber: grid -> raw coords -> one streamed pass over the planes
    %----------------------------------------------------------------------
    Nfibers           = numel(all_results);
    intensity_results = cell(1, Nfibers);
    Dfield            = [];

    for ifiber = 1:Nfibers
        res       = all_results{ifiber};
        depths_um = res.depths_um;
        radius    = res.radius_vox;
        Ndepths   = numel(depths_um);

        % same basis and extent as extractAnnotationCircularSlice
        n   = res.normal_vox(:)' / norm(res.normal_vox);
        tmp = [1 0 0];
        if abs(dot(tmp, n)) > 0.9; tmp = [0 1 0]; end
        e1  = cross(n, tmp);  e1 = e1 / norm(e1);
        e2  = cross(n, e1);   e2 = e2 / norm(e2);

        rvec  = grinSliceGrid(radius, step_vox);
        Ngrid = numel(rvec);
        [G1, G2] = meshgrid(rvec, rvec);

        pts = zeros(Ngrid^2 * Ndepths, 3);         % [AP DV ML], all depths stacked
        for id = 1:Ndepths
            c  = res.center_vox(:)' + depths_um(id) / 10 * n;
            ii = (id - 1) * Ngrid^2 + (1:Ngrid^2);
            pts(ii, :) = c + G1(:) * e1 + G2(:) * e2;
        end

        tic;
        [raw_xyz, Dfield] = atlasPointsToSample(pts(:, [2 1 3]), trstruct, ...
            'UseBspline', opt.UseBspline, 'Dfield', Dfield);
        clear pts;

        vals = sampleRawStack(raw_xyz, tfiles, opts, chans, [Nx Ny Nz]);
        fprintf('  Fiber %d: %d depths × %d² px, %.1f s\n', ...
            ifiber, Ndepths, Ngrid, toc);

        % annotation on the very same grid (nearest neighbour), so the region
        % boundaries line up pixel for pixel with the intensity
        slices_av = {};
        if ~isempty(opt.Annotation)
            slices_av = cell(1, Ndepths);
            for id = 1:Ndepths
                slices_av{id} = extractAnnotationCircularSlice(opt.Annotation, ...
                    res.center_vox(:)', n, radius, depths_um(id) / 10, step_vox);
            end
        end

        % reshape back into per-depth images and disc statistics
        in_circle        = (G1.^2 + G2.^2) <= radius^2;
        slices_int       = cell(1, Ndepths);
        median_intensity = nan(Ndepths, numel(chans));
        mean_intensity   = nan(Ndepths, numel(chans));
        for id = 1:Ndepths
            ii = (id - 1) * Ngrid^2 + (1:Ngrid^2);
            sl = reshape(vals(ii, :), Ngrid, Ngrid, numel(chans));
            slices_int{id} = sl;
            for ichan = 1:numel(chans)
                v = sl(:, :, ichan);
                median_intensity(id, ichan) = median(v(in_circle), 'omitnan');
                mean_intensity(id, ichan)   = mean(v(in_circle), 'omitnan');
            end
        end

        r.slices_int       = slices_int;
        r.median_intensity = median_intensity;
        r.mean_intensity   = mean_intensity;
        r.rvec_arr         = repmat({rvec}, 1, Ndepths);
        r.depths_um        = depths_um;
        r.radius_vox       = radius;
        r.pxsize_um        = pxsize_um;
        r.slices_av        = slices_av;
        intensity_results{ifiber} = r;
    end
end

%==========================================================================

function vals = sampleRawStack(raw_xyz, tfiles, opts, chans, sz)
% Trilinear sample of the raw stack at raw_xyz ([x y z] = [col row plane],
% 1-based), reading each needed plane once and only the crop it needs.
% Points outside the stack get 0, like the registered-volume sampler.
    Nx = sz(1); Ny = sz(2); Nz = sz(3);
    Np = size(raw_xyz, 1);
    vals = zeros(Np, numel(chans), 'single');

    x = raw_xyz(:, 1); y = raw_xyz(:, 2); z = raw_xyz(:, 3);
    valid = x >= 1 & x <= Nx & y >= 1 & y <= Ny & z >= 1 & z <= Nz;
    if ~any(valid); return; end

    % plane below each point (capped so the plane above exists) and z weight
    vid = find(valid);
    z0  = min(floor(z(vid)), max(Nz - 1, 1));
    fz  = single(z(vid) - z0);

    % bucket points by z0
    zlo  = min(z0); zhi = max(z0);
    grp  = z0 - zlo + 1;
    [~, ord] = sort(grp);
    cnt   = accumarray(grp, 1, [zhi - zlo + 1, 1]);
    stops = cumsum(cnt); starts = stops - cnt + 1;
    bucket = @(g) ord(starts(g):stops(g));

    tifftype = getOr(opts, 'tifftype', 'planeperfile');
    multitif = getOr(opts, 'multitiffs', false);

    for ic = 1:numel(chans)
        switch tifftype
            case 'planeperfile'
                readcrop = @(k, rr, cc) planeCrop(imread(tfiles{k}, ...
                    'PixelRegion', {rr, cc}));
            case 'channelperfile'
                if multitif
                    bf = BioformatsImage(tfiles{chans(ic)});
                    readcrop = @(k, rr, cc) bf.getPlane(k, 1, 1);
                else
                    bf = BioformatsImage(tfiles{1});
                    readcrop = @(k, rr, cc) bf.getPlane(k, chans(ic), 1);
                end
            otherwise
                error('extractGRINRawIntensitySlices: unknown tifftype ''%s''.', tifftype);
        end
        fullplane = strcmp(tifftype, 'channelperfile');

        acc = zeros(numel(vid), 1, 'single');
        nplanes = zhi - zlo + 2;
        for k = zlo:min(zhi + 1, Nz)
            g = k - zlo + 1;                    % points with z0 == k use (1 - fz)
            ib = []; ia = [];
            if g <= numel(cnt); ib = bucket(g);     end
            if g >= 2;          ia = bucket(g - 1); end   % z0 == k-1 use fz
            ids = [ib; ia];
            if isempty(ids); continue; end

            xs = x(vid(ids)); ys = y(vid(ids));
            c1 = max(1, floor(min(xs)) - 1); c2 = min(Nx, ceil(max(xs)) + 1);
            r1 = max(1, floor(min(ys)) - 1); r2 = min(Ny, ceil(max(ys)) + 1);

            im = readcrop(k, [r1 r2], [c1 c2]);
            if fullplane; im = im(r1:r2, c1:c2, 1); end
            v  = interp2(single(im), xs - c1 + 1, ys - r1 + 1, 'linear', single(0));

            w  = [1 - fz(ib); fz(ia)];
            acc(ids) = acc(ids) + w .* v;

            if mod(k - zlo, 200) == 0
                fprintf('    %s: plane %d/%d\n', ...
                    sprintf('chan %d', chans(ic)), k - zlo + 1, nplanes);
            end
        end
        vals(vid, ic) = acc;
    end
end

function im = planeCrop(im)
% planeperfile TIFFs are single-channel; drop any extra samples.
    im = im(:, :, 1);
end
