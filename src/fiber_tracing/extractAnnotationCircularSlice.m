function [slice_av, rvec] = extractAnnotationCircularSlice(av, center, normal, radius_vox, depth_offset_vox)
% EXTRACTANNOTATIONCIRCULARSLICE  Sample a circular cross-section from a 3-D
%   annotation volume (nearest-neighbour, suitable for integer atlas IDs).
%
%   [SLICE_AV, RVEC] = extractAnnotationCircularSlice( ...
%       AV, CENTER, NORMAL, RADIUS_VOX, DEPTH_OFFSET_VOX)
%
%   AV               – 3-D annotation volume (integer parcellation IDs),
%                      indexed as AV(AP, DV, ML) in Allen CCF convention.
%   CENTER           – 1x3 circle centre in atlas voxel coords [AP DV ML].
%   NORMAL           – 1x3 unit vector (fibre axis direction).
%   RADIUS_VOX       – radius of the extraction disk in atlas voxels.
%   DEPTH_OFFSET_VOX – depth along NORMAL from CENTER (0 = bottom face,
%                      positive = deeper into the brain).
%
%   SLICE_AV – square 2-D array of parcellation IDs sampled on a regular
%              grid covering the disk (plus a small margin).
%   RVEC     – coordinate vector [vox] for both image axes, centred on 0.
%
%   The sampling resolution is 1 atlas voxel (= 10 µm) by default, with
%   a 5-voxel margin added beyond RADIUS_VOX on each side.

    % Shift centre along fibre axis to the requested depth
    center_d = center + depth_offset_vox * normal;

    %------------------------------------------------------------------
    % Build orthonormal basis for the sampling plane  ⊥  normal
    %------------------------------------------------------------------
    tmp = [1 0 0];
    if abs(dot(tmp, normal)) > 0.9
        tmp = [0 1 0];
    end
    e1 = cross(normal, tmp);   e1 = e1 / norm(e1);   % 1x3
    e2 = cross(normal, e1);    e2 = e2 / norm(e2);   % 1x3

    %------------------------------------------------------------------
    % Sampling grid in the plane
    %------------------------------------------------------------------
    pad   = 5;                                          % extra voxels beyond radius
    rspan = radius_vox + pad;
    Ngrid = 2 * ceil(rspan) + 1;
    rvec  = linspace(-rspan, rspan, Ngrid);

    [G1, G2] = meshgrid(rvec, rvec);                   % each Ngrid × Ngrid

    % 3-D world coordinates of every grid point: shape = (Ngrid^2) × 3
    pts3d = center_d + G1(:) * e1 + G2(:) * e2;       % broadcast: e1/e2 are 1×3

    %------------------------------------------------------------------
    % Nearest-neighbour look-up in the annotation volume
    %------------------------------------------------------------------
    [Sap, Sdv, Sml] = size(av);

    ap = round(pts3d(:, 1));
    dv = round(pts3d(:, 2));
    ml = round(pts3d(:, 3));

    % Clamp to volume bounds (out-of-brain samples → index 1, which is
    % annotation ID 0 / outside-brain by Allen CCF convention)
    ap = max(1, min(Sap, ap));
    dv = max(1, min(Sdv, dv));
    ml = max(1, min(Sml, ml));

    idx      = sub2ind([Sap, Sdv, Sml], ap, dv, ml);
    slice_av = reshape(av(idx), Ngrid, Ngrid);
end
