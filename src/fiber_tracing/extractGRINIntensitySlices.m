function [intensity_results, channames] = extractGRINIntensitySlices(savepath, all_results)
% EXTRACTGRININTENSITYSLICES  Sample fluorescence intensity from registered
%   channel volumes at the circular cross-sections of each GRIN fiber.
%
%   [INTENSITY_RESULTS, CHANNAMES] = extractGRINIntensitySlices(SAVEPATH, ALL_RESULTS)
%
%   Rather than loading full volumes, this function:
%     1. Pre-computes every query point (all fibers × all depths).
%     2. Finds the tight 3-D bounding box of those points.
%     3. Reads only the ML pages that fall inside that box, cropped in
%        AP and DV to the same box, from each  chanNN_*.tif(f)  file.
%     4. Runs interp3 on the resulting small sub-volume.
%
%   For a 500 µm-diameter fiber with 5 depth steps the sub-volume is
%   typically ~100³ voxels vs ~1320×800×1140 for the full atlas volume.
%
%   Inputs
%   ------
%   SAVEPATH     – LightSuite folder that contains  volume_registered/  and
%                  was previously processed by annotateGRINLens.
%   ALL_RESULTS  – 1×Nfibers cell array of fiber result structs loaded from
%                  grin_fiber*_atlas.mat  (fields: center_vox, normal_vox,
%                  radius_vox, depths_um, rvec_arr).
%
%   Outputs
%   -------
%   INTENSITY_RESULTS – 1×Nfibers cell array; each element is a struct:
%       .slices_int   – 1×Ndepths cell array; each is [Ngrid×Ngrid×Nchan]
%                       single-precision intensity array.
%       .rvec_arr     – 1×Ndepths cell array of rvec coordinate vectors.
%       .depths_um    – 1×Ndepths depth vector in µm.
%       .radius_vox   – fiber radius in 10 µm atlas voxels.
%   CHANNAMES  – 1×Nchan cell array of channel name strings.

    %----------------------------------------------------------------------
    % 1. Locate volume_registered folder and channel files
    %----------------------------------------------------------------------
    voldir = fullfile(savepath, 'volume_registered');
    if ~exist(voldir, 'dir')
        error('extractGRINIntensitySlices: volume_registered/ not found in:\n  %s', savepath);
    end

    flist = [dir(fullfile(voldir, 'chan*.tif')); dir(fullfile(voldir, 'chan*.tiff'))];
    if isempty(flist)
        error('extractGRINIntensitySlices: no chan*.tif(f) files found in %s', voldir);
    end

    % Sort by embedded channel number
    channums = zeros(1, numel(flist));
    for k = 1:numel(flist)
        tok = regexp(flist(k).name, '^chan(\d+)', 'tokens', 'once');
        if ~isempty(tok); channums(k) = str2double(tok{1}); end
    end
    [~, sidx] = sort(channums);
    flist     = flist(sidx);
    Nchannels = numel(flist);

    % Parse channel names from filenames  (chanNN_NAME.tif[f] → NAME)
    channames = cell(1, Nchannels);
    for k = 1:Nchannels
        tok = regexp(flist(k).name, '^chan\d+_(.+)\.tiff?$', 'tokens', 'once');
        channames{k} = getOr(tok, 1, sprintf('ch%d', k));
    end

    % Volume shape from first channel
    info_ch1 = imfinfo(fullfile(flist(1).folder, flist(1).name));
    Sap      = info_ch1(1).Height;   % rows  = AP
    Sdv      = info_ch1(1).Width;    % cols  = DV
    Sml      = numel(info_ch1);      % pages = ML

    %----------------------------------------------------------------------
    % 2. Read AP offset  (atlasaplims in regopts.mat or transform_params.mat)
    %    registered-volume AP index = full-CCF AP − ap_offset
    %----------------------------------------------------------------------
    ap_offset    = 0;
    regopts_file = fullfile(savepath, 'regopts.mat');
    tparams_file = fullfile(savepath, 'transform_params.mat');

    if exist(regopts_file, 'file')
        ro = load(regopts_file, 'opts');
        if isfield(ro, 'opts') && isfield(ro.opts, 'atlasaplims')
            ap_offset = ro.opts.atlasaplims(1) - 1;
        end
    end
    if ap_offset == 0 && exist(tparams_file, 'file')
        tp = load(tparams_file, 'atlasaplims');
        if isfield(tp, 'atlasaplims')
            ap_offset = tp.atlasaplims(1) - 1;
        end
    end
    fprintf('  AP offset (atlasaplims - 1) = %d\n', ap_offset);

    %----------------------------------------------------------------------
    % 3. Pre-compute ALL query points (all fibers × all depths)
    %    Store them per (fiber, depth) so we can reassemble results later.
    %----------------------------------------------------------------------
    Nfibers = numel(all_results);

    % Cell array: pts_store{ifiber}{id} = [N×3] in registered-volume coords
    pts_store = cell(1, Nfibers);
    basis_store = cell(1, Nfibers);   % {e1, e2, n} per fiber

    all_ap = [];   % accumulate for global bounding-box
    all_dv = [];
    all_ml = [];

    for ifiber = 1:Nfibers
        center    = all_results{ifiber}.center_vox;
        normal    = all_results{ifiber}.normal_vox;
        depths_um = all_results{ifiber}.depths_um;
        rvec_arr  = all_results{ifiber}.rvec_arr;
        depths_vox = depths_um / 10;
        Ndepths   = numel(depths_um);

        % Orthonormal basis for the cross-section plane
        n   = normal(:)' / norm(normal);
        tmp = [1 0 0];
        if abs(dot(tmp, n)) > 0.9; tmp = [0 1 0]; end
        e1  = cross(n, tmp);  e1 = e1 / norm(e1);
        e2  = cross(n, e1);   e2 = e2 / norm(e2);
        basis_store{ifiber} = struct('n', n, 'e1', e1, 'e2', e2);

        pts_store{ifiber} = cell(1, Ndepths);
        for id = 1:Ndepths
            center_d = center + depths_vox(id) * n;
            rvec     = rvec_arr{id};
            Ngrid    = numel(rvec);
            [G1, G2] = meshgrid(rvec, rvec);
            pts3d    = center_d + G1(:) * e1 + G2(:) * e2;  % (Ngrid²) × 3

            % Convert to registered-volume coordinates
            ap_q = pts3d(:, 1) - ap_offset;
            dv_q = pts3d(:, 2);
            ml_q = pts3d(:, 3);

            pts_store{ifiber}{id} = [ap_q, dv_q, ml_q];

            all_ap = [all_ap; ap_q]; %#ok<AGROW>
            all_dv = [all_dv; dv_q]; %#ok<AGROW>
            all_ml = [all_ml; ml_q]; %#ok<AGROW>
        end
    end

    %----------------------------------------------------------------------
    % 4. Tight bounding box (+1 voxel margin for interpolation)
    %----------------------------------------------------------------------
    ap_min = max(1,    floor(min(all_ap)) - 1);
    ap_max = min(Sap,  ceil(max(all_ap))  + 1);
    dv_min = max(1,    floor(min(all_dv)) - 1);
    dv_max = min(Sdv,  ceil(max(all_dv))  + 1);
    ml_min = max(1,    floor(min(all_ml)) - 1);
    ml_max = min(Sml,  ceil(max(all_ml))  + 1);

    ap_sub = ap_min:ap_max;
    dv_sub = dv_min:dv_max;
    ml_sub = ml_min:ml_max;

    fprintf('  Sub-volume: AP [%d–%d], DV [%d–%d], ML [%d–%d]  →  %d × %d × %d voxels per channel\n', ...
        ap_min, ap_max, dv_min, dv_max, ml_min, ml_max, ...
        numel(ap_sub), numel(dv_sub), numel(ml_sub));

    %----------------------------------------------------------------------
    % 5. Load sub-volumes: only the needed ML pages, cropped in AP and DV
    %    TIFF layout: page k = ML slice k, rows = AP, cols = DV
    %----------------------------------------------------------------------
    fprintf('Reading sub-volumes from disk:\n');
    subvols = cell(1, Nchannels);

    for ichan = 1:Nchannels
        fpath  = fullfile(flist(ichan).folder, flist(ichan).name);
        info   = imfinfo(fpath);
        Sml_f  = numel(info);                  % actual page count in file
        ml_use = ml_sub(ml_sub <= Sml_f);      % clamp to file bounds

        sub = zeros(numel(ap_sub), numel(dv_sub), numel(ml_sub), 'single');
        for ki = 1:numel(ml_use)
            page       = ml_use(ki);
            slice2d    = imread(fpath, 'Index', page, 'Info', info);
            sub(:,:,ki) = single(slice2d(ap_sub, dv_sub));
        end

        subvols{ichan} = sub;
        fprintf('  %s  (%d pages read out of %d)\n', flist(ichan).name, numel(ml_use), Sml_f);
    end

    %----------------------------------------------------------------------
    % 6. Sample sub-volumes with interp3 and reassemble per-fiber results
    %    Sub-volume coords: ap_local = ap_q − ap_min + 1  (etc.)
    %    interp3(V, Xq, Yq, Zq): Xq→cols(DV), Yq→rows(AP), Zq→pages(ML)
    %----------------------------------------------------------------------
    intensity_results = cell(1, Nfibers);

    for ifiber = 1:Nfibers
        depths_um  = all_results{ifiber}.depths_um;
        rvec_arr   = all_results{ifiber}.rvec_arr;
        radius     = all_results{ifiber}.radius_vox;
        Ndepths    = numel(depths_um);

        slices_int       = cell(1, Ndepths);
        median_intensity = nan(Ndepths, Nchannels);  % [Ndepths × Nchannels]

        for id = 1:Ndepths
            pts    = pts_store{ifiber}{id};   % (Ngrid²) × 3  in volume coords
            rvec   = rvec_arr{id};
            Ngrid  = numel(rvec);

            % Translate to sub-volume (1-based) coordinates
            ap_loc = pts(:,1) - (ap_min - 1);
            dv_loc = pts(:,2) - (dv_min - 1);
            ml_loc = pts(:,3) - (ml_min - 1);

            % Clamp to sub-volume bounds
            ap_loc = max(1, min(numel(ap_sub), ap_loc));
            dv_loc = max(1, min(numel(dv_sub), dv_loc));
            ml_loc = max(1, min(numel(ml_sub), ml_loc));

            slice_int = zeros(Ngrid, Ngrid, Nchannels, 'single');
            for ichan = 1:Nchannels
                vals = interp3(subvols{ichan}, dv_loc, ap_loc, ml_loc, 'linear', single(0));
                slice_int(:,:,ichan) = reshape(vals, Ngrid, Ngrid);
            end
            slices_int{id} = slice_int;

            % Median inside the fiber circle (pixels where r ≤ radius_vox)
            [G1, G2]  = meshgrid(rvec, rvec);
            in_circle = (G1.^2 + G2.^2) <= radius^2;
            for ichan = 1:Nchannels
                v = slice_int(:,:,ichan);
                median_intensity(id, ichan) = median(v(in_circle), 'omitnan');
            end
        end

        res_int.slices_int       = slices_int;
        res_int.median_intensity = median_intensity;   % [Ndepths × Nchannels]
        res_int.rvec_arr         = rvec_arr;
        res_int.depths_um        = depths_um;
        res_int.radius_vox       = radius;
        intensity_results{ifiber} = res_int;
    end
end
