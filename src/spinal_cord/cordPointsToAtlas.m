function finalpts = cordPointsToAtlas(inputpts, trstruct, varargin)
%CORDPOINTSTOATLAS Map points from raw spinal cord pixels into atlas voxels.
%
%   FINALPTS = CORDPOINTSTOATLAS(INPUTPTS, TRSTRUCT) transforms the N x M array
%   INPUTPTS, whose first three columns are [x y z] in the *original* sample
%   (the coordinates cell detection produces), into spinal cord atlas voxels.
%   Columns 4 and beyond are descriptors and are carried through untouched.
%   TRSTRUCT is the transform_params struct written by MULTIOBJCORDREGISTRATION.
%
%   A cord goes through more steps than a brain, and all of them have to be
%   undone in order, because most of them threw information away:
%
%     1. Downsampling to the registration resolution. Uses the same convention
%        as imresize, so a point lands where its own pixel landed.
%     2. Permutation of the array dimensions, which put the long axis last.
%     3. Cropping - both the bounding box around the cord and the slices past
%        the point where the brain starts. Points outside the kept range still
%        get transformed, but they are extrapolated and should be treated with
%        suspicion; SANITIZECELLCOORDS drops whatever leaves the annotation.
%     4. Straightening. Every slice was translated and rotated on its own, so
%        this is a per-slice rigid transform, interpolated at the point's own
%        fractional z (CORDSTRAIGHTENPOINTS).
%     5. The B-spline, applied as a displacement field exactly as for a brain.
%     6. The affine that takes the straightened sample onto the atlas.
%     7. The rostrocaudal flip. Registration runs against a flipped atlas when
%        the sample was imaged caudorostral; this puts the result back into the
%        atlas as it is distributed.
%     8. Resampling from the registration grid onto the native atlas grid.
%
%   FINALPTS = CORDPOINTSTOATLAS(..., 'space', 'registration') stops after step
%   7 and returns coordinates in the isotropic registration grid instead of the
%   native atlas grid - the grid the registered volumes are computed on.
%
%   FINALPTS = CORDPOINTSTOATLAS(..., 'Dfield', D) reuses a displacement field
%   already read from disk, which matters when several point sets are
%   transformed in a row. 'skipbspline' (default false) leaves out step 5
%   altogether, for a quick affine-only check.
%
%   See also TRANSFORMCORDPOINTSTOATLAS, MULTIOBJCORDREGISTRATION,
%   CORDSTRAIGHTENPOINTS, GENERATEREGISTEREDCORDVOLUME.

%==========================================================================
p = inputParser;
addParameter(p, 'space', 'atlas', @(x) ischar(x) || isstring(x));
addParameter(p, 'Dfield', [], @isnumeric);
addParameter(p, 'skipbspline', false, @(x) islogical(x) || isscalar(x));
parse(p, varargin{:});
params = p.Results;

space = validatestring(params.space, {'atlas', 'registration'});
%==========================================================================
if isempty(inputpts)
    finalpts = inputpts;
    return
end

registres = trstruct.registrationres(1);       % isotropic, um
pxsize    = trstruct.ori_pxsize;               % [x y z] um in the raw sample
%==========================================================================
% 1. raw sample pixels -> registration volume voxels
% preprocessLightSheetVolume resizes xy by pxsize(1)/registres and z by
% pxsize(3)/registres; imresize maps pixel centres as s*(u - 0.5) + 0.5.
sfac = [pxsize(1) pxsize(1) pxsize(3)] / registres;
v    = sfac .* (double(inputpts(:, 1:3)) - 0.5) + 0.5;
%==========================================================================
% 2. the permutation that put the long axis last, in array dimensions
a = v(:, [2 1 3]);                     % [dim1 dim2 dim3] = [row col plane]
a = a(:, trstruct.how_to_perm);
q = a(:, [2 1 3]);                     % back to [x y z] of the permuted volume
%==========================================================================
% 3. the crop
q = q - [trstruct.samp_ikeepx(1) trstruct.samp_ikeepy(1) trstruct.samp_ikeeplong(1)] + 1;
%==========================================================================
% 4. the per-slice straightening
q = cordStraightenPoints(q, trstruct.straighten);
%==========================================================================
% 5. the B-spline, as a displacement field in registration voxels
if ~params.skipbspline
    regsize_mm = registres * 1e-3;
    Dfield     = params.Dfield;
    if isempty(Dfield)
        Dfield = transformix([], trstruct.tform_bspline_samp20um_to_atlas_20um_px);
        Dfield = permute(Dfield, [2 3 4 1]) / regsize_mm;
    end
    [Sx, Sy, Sz, ~] = size(Dfield);

    dx = interpn(1:Sx, 1:Sy, 1:Sz, Dfield(:,:,:,1), q(:,1), q(:,2), q(:,3), 'linear');
    dy = interpn(1:Sx, 1:Sy, 1:Sz, Dfield(:,:,:,2), q(:,1), q(:,2), q(:,3), 'linear');
    dz = interpn(1:Sx, 1:Sy, 1:Sz, Dfield(:,:,:,3), q(:,1), q(:,2), q(:,3), 'linear');

    % negative sign is extremely important: the field on file runs the other way
    disp_int = -[dx, dy, dz];
    disp_int(isnan(disp_int)) = 0;
    q = q + disp_int;
end
%==========================================================================
% 6. the affine onto the atlas
r = trstruct.tform_affine_samp20um_to_atlas_20um_px.transformPointsForward(q);
%==========================================================================
% 7. undo the rostrocaudal flip of the atlas
if trstruct.tofliprc
    r(:, 3) = trstruct.atlassize(3) + 1 - r(:, 3);
end
%==========================================================================
% 8. onto the native atlas grid
if strcmp(space, 'atlas')
    % the registration atlas was made with imresize3(av, 'Scale', sc), per array
    % dimension, so the inverse mapping is (u - 0.5)/sc + 0.5
    sc  = trstruct.atlasres ./ (registres * [1 1 1]);
    scx = sc([2 1 3]);                 % reorder to [x y z]
    r   = (r - 0.5) ./ scx + 0.5;
end
%==========================================================================
finalpts = [r, double(inputpts(:, 4:end))];
%==========================================================================
end
