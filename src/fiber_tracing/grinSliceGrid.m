function rvec = grinSliceGrid(radius_vox, step_vox)
% GRINSLICEGRID  In-plane grid of the GRIN cross-sections.
%
%   RVEC = grinSliceGrid(RADIUS_VOX, STEP_VOX)
%
%   Single source for the extent of every GRIN cross-section (atlas
%   annotation and intensity, at any resolution): a square of side
%   6 × RADIUS_VOX centred on the axis, i.e. ±3 radii. RVEC is the
%   coordinate vector for both image axes, in 10 µm atlas voxels, spaced
%   STEP_VOX apart (default 1 = 10 µm) and symmetric about 0.

    if nargin < 2 || isempty(step_vox); step_vox = 1; end
    halfwidth = 3 * radius_vox;
    nhalf     = ceil(halfwidth / step_vox - 1e-9);
    rvec      = (-nhalf:nhalf) * step_vox;
end
