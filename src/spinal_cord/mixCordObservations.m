function [idx, val, keepauto] = mixCordObservations(N, useridx, userval, autoval, radius)
%MIXCORDOBSERVATIONS User points, plus predictions from outside their reach.
%
%   [IDX, VAL] = MIXCORDOBSERVATIONS(N, USERIDX, USERVAL, AUTOVAL, RADIUS)
%   assembles the observations that SPINAL_CORD_ALIGNER fits a single smoothing
%   spline to, over a cord of N slices. The user's points (USERVAL at slices
%   USERIDX) are always observations. The automatic prediction AUTOVAL, which
%   covers every slice, is an observation only where the user has not spoken.
%
%   The user speaks for two kinds of slice:
%     * everything within RADIUS of one of their points, and
%     * everything *between* two of their points, however far apart those are.
%
%   The second rule matters as much as the first. Without it two points a few
%   hundred slices apart leave a gap in the middle that the prediction fills, and
%   the fit lurches out to meet it and back - a spike of tens of degrees sitting
%   between two clicks that agree with each other to within one. Once a stretch
%   of cord has been annotated at both ends, it belongs to the user, and RADIUS
%   only decides how far past the outermost points that ownership extends.
%
%   The prediction is an observation on nearly every slice, so a click competing
%   against it directly would barely move the curve; taking it away first is
%   what lets a handful of clicks take over a stretch of cord.
%
%   Pass AUTOVAL empty to fit the user's points alone.
%
%   [..., KEEPAUTO] = MIXCORDOBSERVATIONS(...) also returns the N x 1 logical
%   mask of slices whose prediction survived, which the GUI draws so the
%   stretches the user has taken over are visible as gaps.
%
%   The returned observations are sorted by slice and hold no duplicates: a
%   slice the user clicked is always suppressed for the prediction.
%
%   See also SPINAL_CORD_ALIGNER, FITCORDSERIES, SMOOTHCORDSERIES.

useridx = useridx(:);
userval = userval(:);

keepauto = true(N, 1);
for k = 1:numel(useridx)
    lo = max(1, useridx(k) - radius);
    hi = min(N, useridx(k) + radius);
    keepauto(lo:hi) = false;
end

% everything the user has bracketed is theirs, whatever the gap
if ~isempty(useridx)
    keepauto(min(useridx):max(useridx)) = false;
end

if isempty(autoval)
    keepauto(:) = false;
    idx = useridx;
    val = userval;
    return
end

autoidx = find(keepauto);
idx     = [useridx; autoidx];
val     = [userval;  autoval(autoidx)];

[idx, isort] = sort(idx);
val          = val(isort);

end
