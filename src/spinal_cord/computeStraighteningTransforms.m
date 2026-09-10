function tforms = computeStraighteningTransforms(straighten, target_center, target_orientation_deg)
%COMPUTESTRAIGHTENINGTRANSFORMS Per-slice rigid transforms that straighten the cord.
%
%   TFORMS = COMPUTESTRAIGHTENINGTRANSFORMS(ST) returns one rigidtform2d per
%   slice, which is what TRANSFORMCORDIMAGESLICES applies to warp the sample
%   into the straightened frame. ST is the struct from CORDSTRAIGHTENPARAMS.
%
%   TFORMS = COMPUTESTRAIGHTENINGTRANSFORMS(ALIGN_OUT, TARGET_CENTER, DEG) is
%   also accepted, and builds ST from the raw aligner output first.
%
%   Each transform is the composition
%       1. translate the centre of the section to the origin,
%       2. rotate so the dorsoventral axis reaches the target orientation,
%       3. translate the origin to TARGET_CENTER.
%
%   CORDSTRAIGHTENPOINTS applies exactly this composition to points, evaluated
%   at fractional slice indices; the two must stay in step, so any change here
%   belongs there as well.
%
%   See also CORDSTRAIGHTENPARAMS, CORDSTRAIGHTENPOINTS,
%   TRANSFORMCORDIMAGESLICES, SPINAL_CORD_ALIGNER.

%--------------------------------------------------------------------------
if ~isfield(straighten, 'cx')
    if nargin < 3, target_orientation_deg = 90; end
    straighten = cordStraightenParams(straighten, target_center, target_orientation_deg);
end
%--------------------------------------------------------------------------
Nslices      = straighten.Nslices;
theta_target = deg2rad(straighten.target_orientation_deg);
tx           = straighten.target_center(1);
ty           = straighten.target_center(2);

tforms(Nslices, 1) = rigidtform2d;

for ii = 1:Nslices
    cx = straighten.cx(ii);
    cy = straighten.cy(ii);

    % the fitted angle points towards anterior; the posterior direction is what
    % gets rotated onto the target orientation
    theta_curr = straighten.theta(ii) + pi;
    theta_rot  = wrapToPiLocal(theta_target - theta_curr);

    c = cos(theta_rot);
    s = sin(theta_rot);

    T_to_origin = [1 0 -cx; 0 1 -cy; 0 0 1];
    R           = [c -s 0;  s  c 0;  0 0 1];
    T_to_target = [1 0  tx; 0 1  ty; 0 0 1];

    tforms(ii) = rigidtform2d(T_to_target * R * T_to_origin);
end
%--------------------------------------------------------------------------
end
