function [center, normal] = fitFiberInAtlas(pts, radius_vox)
% FITFIBERSINATLAS  Fit a circle to 3-D atlas-space edge points.
%
%   [CENTER, NORMAL] = fitFiberInAtlas(PTS, RADIUS_VOX)
%
%   PTS        – Nx3 array of 3-D points in Allen CCF 10 µm voxel coords
%                [AP, DV, ML].  These are the fiber-edge points collected
%                across multiple slices by annotateGRINLens.
%   RADIUS_VOX – expected fiber radius in 10 µm atlas voxels (informational;
%                the circle centre is always fitted by least squares).
%
%   CENTER – 1x3 fitted circle centre in atlas voxel coords [AP DV ML].
%   NORMAL – 1x3 unit vector of the fibre axis (normal to the circle
%            plane), oriented so it points along the direction of
%            decreasing AP index (i.e., toward the brain surface in
%            the typical dorsal-implant configuration).
%
%   Algorithm
%   ---------
%   1. PCA on the point cloud: the direction of maximum variance is the
%      cylinder long axis (the fibre axis / NORMAL direction).
%   2. All points are projected onto the plane perpendicular to this axis.
%   3. Algebraic least-squares circle fit in the 2-D plane gives the
%      centre.  The circle centre is then back-projected to 3-D.

    if size(pts, 1) < 4
        error('fitFiberInAtlas: need at least 4 points.');
    end

    %--------------------------------------------------------------
    % 1. Find fibre axis via PCA on the centred cloud
    %--------------------------------------------------------------
    centroid = median(pts, 1);
    pts_c    = pts - centroid;       % Nx3 centred

    % Scatter matrix; eigenvectors in ascending eigenvalue order
    [V, ~] = eig(pts_c' * pts_c);
    % Column with largest eigenvalue = direction of maximum variance = axis
    fiber_axis = V(:, end)';         % 1x3 unit vector

    %--------------------------------------------------------------
    % 2. Build an orthonormal basis for the plane  ⊥  fiber_axis
    %--------------------------------------------------------------
    tmp = [1 0 0];
    if abs(dot(tmp, fiber_axis)) > 0.9
        tmp = [0 1 0];
    end
    e1 = cross(fiber_axis, tmp);   e1 = e1 / norm(e1);  % 1x3
    e2 = cross(fiber_axis, e1);    e2 = e2 / norm(e2);  % 1x3

    %--------------------------------------------------------------
    % 3. Project points onto the plane and fit a circle
    %--------------------------------------------------------------
    % 2-D coordinates in the plane
    coords2d = [pts_c * e1', pts_c * e2'];   % Nx2

    % Algebraic least-squares circle fit:
    %   (x-a)^2 + (y-b)^2 = r^2
    %   linearised: 2ax + 2by + (r^2 - a^2 - b^2) = x^2 + y^2
    %   => A * [a; b; D] = rhs,   D = r^2 - a^2 - b^2
    A   = [2*coords2d, ones(size(coords2d,1), 1)];
    rhs = sum(coords2d.^2, 2);

    params    = A \ rhs;
    center_2d = params(1:2)';    % [a, b] in plane coords

    % Guard against degenerate fit (e.g. all points collinear)
    if ~all(isfinite(center_2d))
        warning('fitFiberInAtlas: degenerate fit – using point-cloud centroid as centre.');
        center_2d = [0 0];
        params(3) = radius_vox^2;
    end

    %--------------------------------------------------------------
    % 4. Back-project centre to 3-D
    %--------------------------------------------------------------
    center = centroid + center_2d(1) * e1 + center_2d(2) * e2;

    % Report fitted radius for reference
    fitted_r = sqrt(max(0, params(3) + center_2d(1)^2 + center_2d(2)^2));
    fprintf('  fitFiberInAtlas: fitted radius = %.1f vox (expected %.1f vox)\n', ...
        fitted_r, radius_vox);

    %--------------------------------------------------------------
    % 5. Orient the normal consistently (toward larger AP index = deeper)
    %--------------------------------------------------------------
    normal = fiber_axis;
    % Ensure normal points toward increasing atlas depth (larger AP index)
    % by checking which extreme point has the larger AP coordinate
    proj = pts * normal';
    projmaxpts = pts(proj > quantile(proj, 0.95),:);
    if median(projmaxpts(:,2)) < centroid(2)
        % The "deep" end has smaller AP → flip
        normal = -normal;
    end
    % --- NEW CODE TO ADD: SHIFT CENTER TO BOTTOM-MOST POINT ---
    % Find the maximum projection along the downward-pointing normal
    depth_projections = (pts - center) * normal';
    bottom_offset = quantile(depth_projections,0.99);
    
    % Shift the center coordinate to the bottom face of the lens
    center = center + bottom_offset * normal;
    %----------------------------------------------------------------------
end
