function atlas_pts = registrationPointsToAtlas(pts_disp, trstruct, registres, savepath, volsize)
% REGISTRATIONPOINTSTOATLAS  Map points clicked on the registration display
%   volume to Allen CCF atlas voxel coordinates.
%
%   ATLAS_PTS = registrationPointsToAtlas(PTS_DISP, TRSTRUCT, REGISTRES, ...
%                                         SAVEPATH, VOLSIZE)
%
%   Points are picked on the permuted registration display volume shown by
%   annotateGRINLens / annotateNeuropixelsProbes (i.e. permuteBrainVolume
%   applied to chan_X_sample_register_20um.tif).  To stay consistent with the
%   rest of LightSuite, this function does NOT re-implement the sample→atlas
%   transform: it converts the clicked points back into full-resolution sample
%   coordinates and then defers to transformPointsToAtlas / coreTransform — the
%   same validated similarity→affine→B-spline chain used for cell mapping.
%
%   Inputs
%   ------
%   PTS_DISP  – N×3 voxel coords in the permuted display volume. Columns are
%               [dim1 dim2 dim3], matching gui_data.all_points in the GUIs.
%   TRSTRUCT  – struct loaded from transform_params.mat (needs how_to_perm,
%               ori_pxsize, ori_size, the bspline path and the affine, as
%               required by coreTransform).
%   REGISTRES – registration resolution in µm (regopts.registres, e.g. 20).
%   SAVEPATH  – LightSuite output folder (used to resolve the bspline file).
%   VOLSIZE   – size of the permuted display volume, size(voldisp); needed to
%               invert the display flips.
%
%   Output
%   ------
%   ATLAS_PTS – N×3 atlas points in the same column order coreTransform
%               returns. ATLAS_PTS(:, [2 1 3]) gives [AP DV ML] voxel
%               coordinates suitable for indexing annotation_10.nii.gz.

    htp     = trstruct.how_to_perm;
    absperm = abs(htp);

    %------------------------------------------------------------------
    % 1. Undo the display flips applied by permuteBrainVolume (it flips the
    %    permuted volume along every dim where how_to_perm < 0).
    %------------------------------------------------------------------
    pd = pts_disp;
    for d = 1:3
        if htp(d) < 0
            pd(:, d) = volsize(d) + 1 - pd(:, d);
        end
    end

    %------------------------------------------------------------------
    % 2. Undo the display permute. permuteBrainVolume does
    %    voldisp = permute(volraw, absperm), so display dim d corresponds to
    %    sample (volraw) dim absperm(d).
    %------------------------------------------------------------------
    s20 = zeros(size(pd), 'like', pd);
    s20(:, absperm) = pd;                 % sample-orientation 20 µm voxels [Y X Z]

    %------------------------------------------------------------------
    % 3. Convert 20 µm sample voxels to full-resolution sample coordinates
    %    [x y z] — the input expected by transformPointsToAtlas/coreTransform.
    %    coreTransform multiplies (pts-1) by ori_pxsize, so this reproduces the
    %    clicked point's physical position exactly regardless of ori_pxsize.
    %------------------------------------------------------------------
    s20_xyz    = s20(:, [2 1 3]);                                  % [x y z]
    sample_pts = (s20_xyz - 1) .* (registres ./ trstruct.ori_pxsize) + 1;

    %------------------------------------------------------------------
    % 4. Resolve the bspline file relative to savepath (robust to moved
    %    folders), then defer to the validated transform.
    %------------------------------------------------------------------
    truse = trstruct;
    [~, bn, be] = fileparts(trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    bpath = fullfile(savepath, [bn be]);
    if exist(bpath, 'file')
        truse.tform_bspline_samp20um_to_atlas_20um_px = bpath;
    end

    atlas_pts = transformPointsToAtlas(sample_pts, 'transform_params', truse);
    atlas_pts = atlas_pts(:, 1:3);
end
