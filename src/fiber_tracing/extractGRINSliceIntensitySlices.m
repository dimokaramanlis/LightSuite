function [intensity_results, channames] = extractGRINSliceIntensitySlices(savepath, all_results, varargin)
% EXTRACTGRINSLICEINTENSITYSLICES  Sample the GRIN cross-sections of a
%   slice-module brain from the aligned sections at px_process resolution.
%
%   [INTENSITY_RESULTS, CHANNAMES] = extractGRINSliceIntensitySlices(SAVEPATH, ALL_RESULTS)
%
%   Slice-module counterpart of extractGRINRawIntensitySlices: same planes,
%   extent and outputs, but every grid point is mapped atlas -> section pixel
%   of volume_aligned/ (sliceAtlasPointsToSection, the chain that builds
%   volume_registered/) and bilinearly sampled there. Each point takes the
%   nearest section, as volume_registered/ does; points beyond the first or
%   last section are NaN.
%
%   Name-value options
%   ------------------
%   'PixelSize'  ([])  grid spacing in µm. Default: sliceinfo.px_process.
%   'Channels'   ([])  channels (chanNN files) to sample. Default: all.
%   'Annotation' ([])  annotation_10 volume; when given, the parcellation is
%                      also sampled on the same fine grid (.slices_av).
%
%   Outputs as extractGRINRawIntensitySlices.

    p = inputParser;
    addParameter(p, 'PixelSize',  [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'Channels',   [], @isnumeric);
    addParameter(p, 'Annotation', [], @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    trstruct  = load(fullfile(savepath, 'transform_params.mat'));
    s         = load(fullfile(savepath, 'sliceinfo.mat'), 'sliceinfo');
    sliceinfo = s.sliceinfo;

    [files, channames] = sliceChannelFiles(savepath, sliceinfo);
    chans = opt.Channels;
    if isempty(chans); chans = 1:numel(files); end
    files     = files(chans);
    channames = channames(chans);

    pxsize_um = opt.PixelSize;
    if isempty(pxsize_um); pxsize_um = sliceinfo.px_process; end
    step_vox = pxsize_um / 10;

    fprintf('Sampling aligned sections at %.2f µm (%d channel(s)):\n', ...
        pxsize_um, numel(chans));

    Nfibers           = numel(all_results);
    intensity_results = cell(1, Nfibers);
    Dfields           = {};

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
        [isec, rc, Dfields] = sliceAtlasPointsToSection(pts, trstruct, sliceinfo, savepath, Dfields);
        clear pts;
        vals = sampleSections(isec, rc, files);
        fprintf('  Fiber %d: %d depths × %d² px, sections %s, %.1f s\n', ifiber, ...
            Ndepths, Ngrid, mat2str(unique(isec(~isnan(isec)))'), toc);

        slices_av = {};
        if ~isempty(opt.Annotation)
            slices_av = cell(1, Ndepths);
            for id = 1:Ndepths
                slices_av{id} = extractAnnotationCircularSlice(opt.Annotation, ...
                    res.center_vox(:)', n, radius, depths_um(id) / 10, step_vox);
            end
        end

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

function vals = sampleSections(isec, rc, files)
% bilinear sample of each needed section, reading only the crop it needs;
% NaN outside the sections (no section, or off the section image)
    vals = nan(numel(isec), numel(files), 'single');
    for ic = 1:numel(files)
        info = imfinfo(files{ic});
        H = info(1).Height;  W = info(1).Width;
        for k = unique(isec(~isnan(isec)))'
            ii = find(isec == k);
            r  = rc(ii, 1);  c = rc(ii, 2);
            on = r >= 1 & r <= H & c >= 1 & c <= W;
            if ~any(on); continue; end
            r1 = max(1, floor(min(r(on))) - 1); r2 = min(H, ceil(max(r(on))) + 1);
            c1 = max(1, floor(min(c(on))) - 1); c2 = min(W, ceil(max(c(on))) + 1);
            im = imread(files{ic}, 'Index', k, 'Info', info, 'PixelRegion', {[r1 r2], [c1 c2]});
            vals(ii(on), ic) = interp2(single(im), c(on) - c1 + 1, r(on) - r1 + 1, 'linear');
        end
    end
end

function [files, names] = sliceChannelFiles(savepath, sliceinfo)
% chanNN_*.tif(f) of volume_aligned/ (next to SAVEPATH, else sliceinfo.slicevolfin)
    voldir = fullfile(savepath, 'volume_aligned');
    if ~exist(voldir, 'dir'); voldir = sliceinfo.slicevolfin; end
    flist = [dir(fullfile(voldir, 'chan*.tif')); dir(fullfile(voldir, 'chan*.tiff'))];
    if isempty(flist)
        error('extractGRINSliceIntensitySlices:notFound', ...
            'No chanNN_*.tif(f) section files in %s', voldir);
    end
    nums = cellfun(@(s) str2double(regexp(s, '^chan(\d+)', 'tokens', 'once')), {flist.name});
    [~, o] = sort(nums);
    flist  = flist(o);
    files  = fullfile({flist.folder}, {flist.name});
    names  = cellfun(@(s) regexprep(s, '^chan\d+_|\.tiff?$', ''), {flist.name}, 'UniformOutput', false);
end
