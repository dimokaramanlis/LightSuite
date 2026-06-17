function atlas_pts = registrationPointsToAtlas(pts_disp, trstruct, registres, savepath)
% REGISTRATIONPOINTSTOATLAS  Map points clicked on the registration volume to
%   Allen CCF atlas voxel coordinates.
%
%   ATLAS_PTS = registrationPointsToAtlas(PTS_DISP, TRSTRUCT, REGISTRES, SAVEPATH)
%
%   Transforms points that were picked on the permuted LightSuite registration
%   display volume (the same volume shown by annotateGRINLens and
%   annotateNeuropixelsProbes) into Allen CCF 10 µm atlas voxels, applying the
%   bspline deformation field followed by the affine transform produced by the
%   registration pipeline.  This is the shared coordinate transform reused by
%   both the GRIN lens and the Neuropixels probe tracing tools.
%
%   Inputs
%   ------
%   PTS_DISP  – N×3 voxel coordinates in the permuted registration display
%               volume.  Columns are [slice(AP), row(DV), col(ML)], i.e. the
%               same convention as gui_data.all_points in the annotation GUIs.
%   TRSTRUCT  – struct loaded from transform_params.mat.  Must contain
%               how_to_perm, ori_pxsize,
%               tform_bspline_samp20um_to_atlas_20um_px (a file path) and
%               tform_affine_samp20um_to_atlas_10um_px (an affine transform).
%   REGISTRES – registration resolution in µm (regopts.registres, e.g. 20).
%   SAVEPATH  – LightSuite output folder; used to resolve the bspline
%               deformation file relative to the (possibly moved) folder.
%
%   Output
%   ------
%   ATLAS_PTS – N×3 points in Allen CCF atlas space.  The column order matches
%               the GRIN pipeline: ATLAS_PTS(:, [2 1 3]) gives [AP DV ML] voxel
%               coordinates suitable for indexing annotation_10.nii.gz.

    regsize_mm = registres * 1e-3;

    %------------------------------------------------------------------
    % 1. Convert display-volume voxels to registration (20 µm) voxels.
    %    The ori_pxsize factors cancel; this reduces to a [2 1 3] reorder
    %    because the display volume is already at the registration
    %    resolution.  Kept explicit to mirror the original GRIN pipeline.
    %------------------------------------------------------------------
    absperm = abs(trstruct.how_to_perm);
    facmult = registres ./ trstruct.ori_pxsize;
    facmult = facmult(absperm);
    ptsuse  = pts_disp .* facmult .* trstruct.ori_pxsize(absperm) * 1e-3;
    pts     = ptsuse(:, [2 1 3]) / regsize_mm;

    %------------------------------------------------------------------
    % 2. Apply the bspline deformation field (resolved relative to savepath)
    %------------------------------------------------------------------
    [~, trname, trext] = fileparts(trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    bsplinepath = fullfile(savepath, [trname trext]);
    if ~exist(bsplinepath, 'file')
        bsplinepath = trstruct.tform_bspline_samp20um_to_atlas_20um_px;
    end

    Dfield = transformix([], bsplinepath);
    Dfield = permute(Dfield, [2 3 4 1]) / regsize_mm;
    [Sx, Sy, Sz, ~] = size(Dfield);

    Xgv = 1:Sx; Ygv = 1:Sy; Zgv = 1:Sz;
    dx = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,1), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dy = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,2), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dz = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,3), pts(:,1), pts(:,2), pts(:,3), 'linear');

    interpolated_displacements = -[dx, dy, dz];
    interpolated_displacements(isnan(interpolated_displacements)) = 0;

    %------------------------------------------------------------------
    % 3. Apply the affine transform to land in 10 µm atlas voxels
    %------------------------------------------------------------------
    atlas_pts = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward( ...
        pts + interpolated_displacements);
    atlas_pts = atlas_pts(:, 1:3);
end
