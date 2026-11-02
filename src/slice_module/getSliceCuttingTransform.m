function [tformuse, anglesource] = getSliceCuttingTransform(tformrigid, Ratlas, procpath, useautoangle)
% getSliceCuttingTransform returns the atlas -> sample rigid transform used
% by the slice control point GUI and by slice registration.
%
% The cutting angle estimated by the automated rigid alignment is not
% reliable, so by default it is removed (cutting angle = 0, i.e. atlas
% coronal planes map onto sample slices). The in-plane rotation and the
% position of the atlas centre are kept.
%
% Priority:
%   1. cutting_angle_data.mat with at least one stored slice -> user angle
%   2. useautoangle = true -> automated alignment as is (legacy behaviour)
%   3. otherwise -> zero cutting angle
%
% matchControlPointsInSlices and registerSlicesToAtlas must call this with
% the same inputs, since control points are stored in the space of the
% transformed atlas.
%
% anglesource is 'user', 'auto' or 'zero'.

if nargin < 4 || isempty(useautoangle)
    useautoangle = false;
end

cutting_angle_data = loadCuttingAngleData(procpath);

if ~isempty(cutting_angle_data)
    fprintf('Using the cutting angle from cutting_angle_data.mat\n');
    tformuse    = applyAngleToTransform(tformrigid, Ratlas, cutting_angle_data);
    anglesource = 'user';
elseif useautoangle
    fprintf('Using the cutting angle of the automated alignment\n');
    tformuse    = tformrigid;
    anglesource = 'auto';
else
    fprintf('No cutting angle defined, assuming 0 (run determineCuttingAngleGUI to set it)\n');
    tformuse    = removeCuttingAngle(tformrigid, Ratlas);
    anglesource = 'zero';
end
end

function cutting_angle_data = loadCuttingAngleData(procpath)
% Returns [] if there is no file or it has no stored slices.
cutting_angle_data = [];
angle_data_path    = fullfile(procpath, 'cutting_angle_data.mat');
if ~exist(angle_data_path, 'file')
    return;
end
data = load(angle_data_path);
if ~isfield(data, 'cutting_angle_data') || ...
        ~isfield(data.cutting_angle_data, 'saved_slice_vectors')
    warning('%s has no saved slice vectors, ignoring it.', angle_data_path);
    return;
end
vecs     = data.cutting_angle_data.saved_slice_vectors;
hasvalid = any(cellfun(@(v) numel(v) == 3 && all(~isnan(v)), vecs));
if ~hasvalid
    warning('%s has no stored slices (press Enter in determineCuttingAngleGUI), ignoring it.', ...
        angle_data_path);
    return;
end
cutting_angle_data = data.cutting_angle_data;
end

function tformout = removeCuttingAngle(tformin, Ratlas)
% Rotate the transform so that the atlas AP axis (world Y) maps onto the
% sample slice axis (world Y). This is the smallest rotation that removes
% the tilt (swing-twist decomposition), so the rotation about AP (in-plane)
% is kept. The atlas centre maps to the same sample point as before.

A  = tformin.A;
R  = A(1:3, 1:3);
T  = A(1:3, 4);
ey = [0; 1; 0];

v  = R * ey;                 % where the atlas AP axis ends up
v  = v / norm(v);
target = sign(v(2)) * ey;    % keep the AP direction (no flips)
if target(2) == 0, target = ey; end

k  = cross(v, target);
s  = norm(k);
c  = dot(v, target);
if s > 1e-12
    k = k / s;
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    Q = eye(3) + s*K + (1 - c)*(K*K);
else
    Q = eye(3);
end
Rnew = Q * R;

atlas_center = [mean(Ratlas.XWorldLimits); mean(Ratlas.YWorldLimits); mean(Ratlas.ZWorldLimits)];
Tnew = T + R*atlas_center - Rnew*atlas_center;

Anew = eye(4);
Anew(1:3, 1:3) = Rnew;
Anew(1:3, 4)   = Tnew;
tformout = rigidtform3d(Anew);
end
