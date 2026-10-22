function [out, D] = brainPointTransform(pts, trstruct, direction, varargin)
%BRAINPOINTTRANSFORM The point transform between a brain sample and the atlas.
%
%   OUT = BRAINPOINTTRANSFORM(PTS, TRSTRUCT, DIRECTION) maps the N x M array
%   PTS between the three coordinate systems of a light-sheet brain. Only the
%   first three columns are coordinates; columns 4:end are descriptors and
%   are carried through untouched. TRSTRUCT is the struct saved as
%   transform_params.mat by multiobjRegistration. This is the ONLY place the
%   brain chain is written down; every other point mapping calls it.
%
%   Coordinate systems (all 1-based, all [x y z] = [col row plane]):
%     'raw'    full-resolution sample pixels, the 'cell_locations' convention.
%     'reg'    voxels of the permuted 20 um registration volume, i.e. of
%              permuteBrainVolume(sample_register_20um.tif, how_to_perm).
%     'atlas'  voxels of the 10 um Allen CCFv3 volumes, as [dim2 dim1 dim3]
%              of annotation_10 (= [DV AP ML]); round and swap to index it.
%
%   DIRECTION is one of 'raw2atlas', 'atlas2raw', 'raw2reg', 'reg2raw',
%   'reg2atlas', 'atlas2reg'.
%
%   The chain, raw -> atlas:
%     1. Downsample: v = s .* (u - 0.5) + 0.5 with s = ori_pxsize / 20 um,
%        the pixel-centre convention of imresize/imresize3, which made the
%        registration volume (preprocessLightSheetVolume).
%     2. Permute and flip exactly as permuteBrainVolume, on the 20 um grid
%        (regvolsize), so a flip maps voxel i to N + 1 - i.
%     3. B-spline: the FORWARD file (bspline_atlas_to_samp_20um.txt) moves a
%        sample point onto the affinely warped atlas, q' = q + D(q), with D
%        read at the point's own position. The inverse file stored in
%        transform_params is for images and must not be used for points.
%     4. Affine onto the 10 um atlas, applied to 1-based voxel coordinates,
%        the convention of the imwarp call it was fitted and used with.
%   atlas -> raw runs the same steps backwards; the B-spline step is
%   inverted exactly (fixed-point iteration), so a round trip closes to
%   ~1e-4 voxels.
%
%   Options (name-value):
%     'UseBspline' (true)  false skips step 3 (affine only).
%     'Dfield'     ([])    field from a previous call (second output), so
%                          transformix runs once for many point sets.
%     'savepath'   ('')    registration folder, searched first for the
%                          forward B-spline file.
%     'MaxIter', 'Tol'     passed to the B-spline inversion.
%
%   [OUT, D] = ... also returns the displacement field for reuse.
%
%   See also TRANSFORMPOINTSTOATLAS, ATLASPOINTSTOSAMPLE, READBSPLINEFIELD,
%   APPLYBSPLINEFIELD, PERMUTEBRAINVOLUME.

ip = inputParser;
addParameter(ip, 'UseBspline', true, @(x) islogical(x) || isscalar(x));
addParameter(ip, 'Dfield',     [],   @isnumeric);
addParameter(ip, 'savepath',   '',   @(x) ischar(x) || isstring(x));
addParameter(ip, 'MaxIter',    50,   @(x) isscalar(x) && x >= 0);
addParameter(ip, 'Tol',        1e-4, @(x) isscalar(x) && x > 0);
parse(ip, varargin{:});
opt = ip.Results;

direction = validatestring(direction, ...
    {'raw2atlas', 'atlas2raw', 'raw2reg', 'reg2raw', 'reg2atlas', 'atlas2reg'});

D = opt.Dfield;
if isempty(pts)
    out = pts;
    return;
end
if size(pts, 2) < 3
    error('brainPointTransform:badInput', 'Points need at least 3 columns [x y z].');
end

xyz   = double(pts(:, 1:3));
regum = trstruct.atlasres * 2;                         % registration voxel, um

needfield = opt.UseBspline && contains(direction, 'atlas');
if needfield && isempty(D)
    D = readBsplineField(trstruct, regum * 1e-3, opt.savepath);
end
invopts = {'MaxIter', opt.MaxIter, 'Tol', opt.Tol};

switch direction
    case 'raw2reg'
        c = rawToReg(xyz, trstruct, regum);
    case 'reg2raw'
        c = regToRaw(xyz, trstruct, regum);
    case 'reg2atlas'
        c = regToAtlas(xyz, trstruct, D, opt.UseBspline);
    case 'atlas2reg'
        c = atlasToReg(xyz, trstruct, D, opt.UseBspline, invopts);
    case 'raw2atlas'
        c = regToAtlas(rawToReg(xyz, trstruct, regum), trstruct, D, opt.UseBspline);
    case 'atlas2raw'
        c = regToRaw(atlasToReg(xyz, trstruct, D, opt.UseBspline, invopts), trstruct, regum);
end

out = [c, double(pts(:, 4:end))];
if isfloat(pts)
    out = cast(out, 'like', pts);
end
end

%==========================================================================
function q = rawToReg(u, trstruct, regum)
s = trstruct.ori_pxsize / regum;                       % [x y z]
v = s .* (u - 0.5) + 0.5;                              % unpermuted reg [x y z]
a = v(:, [2 1 3]);                                     % [row col plane]
[perm, flips, sz] = permInfo(trstruct, s);
b = a(:, perm);
for d = find(flips)
    b(:, d) = sz(perm(d)) + 1 - b(:, d);
end
q = b(:, [2 1 3]);
end

%--------------------------------------------------------------------------
function u = regToRaw(q, trstruct, regum)
s = trstruct.ori_pxsize / regum;
[perm, flips, sz] = permInfo(trstruct, s);
b = q(:, [2 1 3]);
for d = find(flips)
    b(:, d) = sz(perm(d)) + 1 - b(:, d);
end
a = zeros(size(b));
a(:, perm) = b;
v = a(:, [2 1 3]);
u = (v - 0.5) ./ s + 0.5;
end

%--------------------------------------------------------------------------
function r = regToAtlas(q, trstruct, D, usebspline)
if usebspline
    q = applyBsplineField(D, q, 'forward');
end
r = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward(q);
end

%--------------------------------------------------------------------------
function q = atlasToReg(r, trstruct, D, usebspline, invopts)
q = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsInverse(r);
if usebspline
    q = applyBsplineField(D, q, 'inverse', invopts{:});
end
end

%--------------------------------------------------------------------------
function [perm, flips, sz] = permInfo(trstruct, s)
% sz: size of the UNPERMUTED registration volume [Ny Nx Nz] at 20 um
perm  = abs(trstruct.how_to_perm);
flips = trstruct.how_to_perm < 0;
if isfield(trstruct, 'regvolsize') && ~isempty(trstruct.regvolsize)
    sz = trstruct.regvolsize(1:3);
else
    % imresize/imresize3 output size for a scale factor
    sz = ceil(trstruct.ori_size .* s([2 1 3]));
end
end
