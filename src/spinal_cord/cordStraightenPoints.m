function ptsout = cordStraightenPoints(pts, straighten, direction)
%CORDSTRAIGHTENPOINTS Straighten (or unstraighten) points along the cord.
%
%   PTSOUT = CORDSTRAIGHTENPOINTS(PTS, ST) applies to the N x 3 points PTS the
%   same per-slice rigid transform that COMPUTESTRAIGHTENINGTRANSFORMS applies
%   to the images. PTS are [x y z] in the cropped registration volume; the third
%   column is untouched, since straightening only moves things within a plane.
%
%   The difference from the image path is that the centre and the angle are
%   interpolated linearly at the point's own (fractional) z, instead of being
%   taken from the nearest slice. A detected cell rarely sits exactly on a slice
%   centre, and rounding it there would shift it by up to half a slice of twist.
%   Points outside the slice range are clamped to the first or last slice.
%
%   PTSOUT = CORDSTRAIGHTENPOINTS(PTS, ST, 'inverse') goes the other way, from
%   the straightened frame back to the cropped registration volume.
%
%   See also CORDSTRAIGHTENPARAMS, COMPUTESTRAIGHTENINGTRANSFORMS,
%   CORDPOINTSTOATLAS.

if nargin < 3 || isempty(direction)
    direction = 'forward';
end

ptsout = pts;
if isempty(pts)
    return
end

N = straighten.Nslices;
z = double(pts(:, 3));
z = min(max(z, 1), N);

if N > 1
    slicevec = (1:N)';
    cx    = interp1(slicevec, straighten.cx,        z, 'linear', 'extrap');
    cy    = interp1(slicevec, straighten.cy,        z, 'linear', 'extrap');
    thraw = interp1(slicevec, straighten.thetaunwr, z, 'linear', 'extrap');
else
    cx    = repmat(straighten.cx(1), size(z));
    cy    = repmat(straighten.cy(1), size(z));
    thraw = repmat(straighten.thetaunwr(1), size(z));
end

theta_target = deg2rad(straighten.target_orientation_deg);
theta_rot    = wrapToPiLocal(theta_target - (thraw + pi));
c            = cos(theta_rot);
s            = sin(theta_rot);

tx = straighten.target_center(1);
ty = straighten.target_center(2);

x = double(pts(:, 1));
y = double(pts(:, 2));

switch lower(direction)
    case 'forward'
        dx = x - cx;
        dy = y - cy;
        ptsout(:, 1) = tx + c .* dx - s .* dy;
        ptsout(:, 2) = ty + s .* dx + c .* dy;
    case 'inverse'
        dx = x - tx;
        dy = y - ty;
        ptsout(:, 1) = cx + c .* dx + s .* dy;
        ptsout(:, 2) = cy - s .* dx + c .* dy;
    otherwise
        error('cordStraightenPoints:badDirection', ...
            'direction must be ''forward'' or ''inverse'', not ''%s''.', direction);
end

end
