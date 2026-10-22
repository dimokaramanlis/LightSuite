function [raw_xyz, Dfield] = atlasPointsToSample(atlas_pts, trstruct, varargin)
%ATLASPOINTSTOSAMPLE Map atlas points back to full-resolution raw sample pixels.
%
%   RAW_XYZ = ATLASPOINTSTOSAMPLE(ATLAS_PTS, TRSTRUCT) is the exact inverse of
%   the cell transform (transformPointsToAtlas): it takes atlas points in the
%   column order that transform RETURNS and gives back raw pixel coordinates
%   [x y z] = [col row plane], 1-based, i.e. the 'cell_locations' convention.
%   Feeding RAW_XYZ to transformPointsToAtlas returns ATLAS_PTS (to ~1e-4
%   registration voxels).
%
%   ATLAS_PTS columns are [dim2 dim1 dim3] of annotation_10 (= [DV AP ML]);
%   a point stored as [AP DV ML] goes in as pts(:, [2 1 3]).
%
%   It is the same as brainPointTransform(ATLAS_PTS, TRSTRUCT, 'atlas2raw').
%
%   Name-value options:
%       'UseBspline' (true)  include the b-spline warp. false = affine only.
%       'Dfield'     ([])    displacement field from a previous call (second
%                            output), so transformix runs only once.
%       'savepath'   ('')    registration folder holding the B-spline files.
%       'MaxIter'    (50)    iterations of the B-spline inversion.
%       'Tol'        (1e-4)  stop when no point moves more than this (20 um px).
%
%   [RAW_XYZ, DFIELD] = ... also returns the displacement field for reuse.
%
%   See also BRAINPOINTTRANSFORM, TRANSFORMPOINTSTOATLAS, RENDERREGIONRAWMIP.

[raw_xyz, Dfield] = brainPointTransform(double(atlas_pts(:, 1:3)), trstruct, ...
    'atlas2raw', varargin{:});
end
