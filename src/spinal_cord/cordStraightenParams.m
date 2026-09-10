function st = cordStraightenParams(align_out, target_center, target_orientation_deg)
%CORDSTRAIGHTENPARAMS Pack the aligner output into a straightening description.
%
%   ST = CORDSTRAIGHTENPARAMS(ALIGN_OUT, TARGET_CENTER) turns the fit saved by
%   SPINAL_CORD_ALIGNER into the compact struct that defines, for every slice,
%   the rigid transform that untwists and unbends the cord: move the centre of
%   the section to TARGET_CENTER and rotate its dorsoventral axis to a fixed
%   orientation.
%
%   ST = CORDSTRAIGHTENPARAMS(..., TARGET_ORIENTATION_DEG) sets that orientation
%   (default 90 degrees, i.e. the posterior direction points down the image, so
%   anterior ends up at the top).
%
%   This struct - and not the array of rigidtform2d objects - is what gets
%   stored in transform_params.mat, because it can be evaluated at a fractional
%   slice index. Volumes are warped slice by slice, but a detected cell sits
%   between slices, and interpolating the centre and the angle is the only way
%   to place it consistently with the image it came from.
%
%   OUTPUT
%     st - struct with fields
%            cx, cy    - N x 1 centre of every section (pixels).
%            theta     - N x 1 dorsoventral angle (rad), pointing from the
%                        posterior towards the anterior side.
%            thetaunwr - N x 1 the same angle unwrapped, so it can be
%                        interpolated across the +/-pi seam.
%            target_center, target_orientation_deg - the straightened frame.
%            Nslices   - number of slices.
%
%   See also SPINAL_CORD_ALIGNER, COMPUTESTRAIGHTENINGTRANSFORMS,
%   CORDSTRAIGHTENPOINTS.

if nargin < 3 || isempty(target_orientation_deg)
    target_orientation_deg = 90;
end

cx    = align_out.fit_x(:);
cy    = align_out.fit_y(:);
theta = align_out.fit_theta(:);

% the aligner leaves NaNs wherever nothing constrained the fit
if any(isnan(cx)) || any(isnan(cy))
    [cx, cy] = deal(fillmissing(cx, 'nearest'), fillmissing(cy, 'nearest'));
end
if any(isnan(theta))
    theta = fillmissing(theta, 'nearest');
end

st = struct();
st.cx        = cx;
st.cy        = cy;
st.theta     = theta;
st.thetaunwr = unwrap(theta);
st.target_center           = reshape(target_center, 1, 2);
st.target_orientation_deg  = target_orientation_deg;
st.Nslices   = numel(cx);

end
