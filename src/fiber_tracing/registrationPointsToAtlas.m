function atlas_pts = registrationPointsToAtlas(pts_disp, trstruct, registres, savepath, volsize) %#ok<INUSD>
% REGISTRATIONPOINTSTOATLAS  Map points clicked on the registration display
%   volume to Allen CCF atlas voxel coordinates.
%
%   ATLAS_PTS = registrationPointsToAtlas(PTS_DISP, TRSTRUCT, REGISTRES, ...
%                                         SAVEPATH, VOLSIZE)
%
%   Points are picked on the permuted registration display volume shown by
%   annotateGRINLens / annotateNeuropixelsProbes (i.e. permuteBrainVolume
%   applied to chan_X_sample_register_20um.tif). For light-sheet brains this
%   is the registration grid itself, so the points go straight through
%   brainPointTransform ('reg2atlas'), the same chain used for cells.
%
%   Inputs
%   ------
%   PTS_DISP  – N×3 voxel coords in the permuted display volume. Columns are
%               [dim1 dim2 dim3], matching gui_data.all_points in the GUIs.
%   TRSTRUCT  – struct loaded from transform_params.mat.
%   REGISTRES – registration resolution in µm (regopts.registres, e.g. 20);
%               used for slice-module brains only.
%   SAVEPATH  – LightSuite output folder (used to find the bspline file).
%   VOLSIZE   – unused; kept so existing callers keep working.
%
%   Output
%   ------
%   ATLAS_PTS – N×3 atlas points in the column order brainPointTransform
%               returns. ATLAS_PTS(:, [2 1 3]) gives [AP DV ML] voxel
%               coordinates suitable for indexing annotation_10.nii.gz.

    %------------------------------------------------------------------
    % Slice-module brains: the display volume is volume_registered/ (atlas
    % space, AP cropped at atlasaplims) block-averaged by ds, with dims
    % [AP DV ML] (see loadTracingRegistration). Display voxel i spans
    % registered voxels (i-1)*ds+1 .. i*ds, centred at (i-0.5)*ds + 0.5.
    %------------------------------------------------------------------
    if isSliceRegistration(trstruct)
        ds  = registres / trstruct.atlasres;
        reg = (pts_disp - 0.5) * ds + 0.5;                 % [AP DV ML], cropped
        reg(:, 1) = reg(:, 1) + trstruct.atlasaplims(1) - 1;
        atlas_pts = reg(:, [2 1 3]);                       % brainPointTransform order
        return;
    end

    %------------------------------------------------------------------
    % Light-sheet brains: the display volume IS the permuted registration
    % volume, so [dim1 dim2 dim3] display voxels are registration voxels and
    % only need reordering to [x y z] for brainPointTransform. VOLSIZE is no
    % longer needed (the flips are already part of that grid).
    %------------------------------------------------------------------
    reg_xyz   = pts_disp(:, [2 1 3]);
    atlas_pts = brainPointTransform(reg_xyz, trstruct, 'reg2atlas', 'savepath', savepath);
    atlas_pts = atlas_pts(:, 1:3);
end
