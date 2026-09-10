function yfit = smoothCordSeries(N, obsidx, obsval, lambda, obswt)
%SMOOTHCORDSERIES Smooth interpolation of sparse observations along the cord.
%
%   YFIT = SMOOTHCORDSERIES(N, OBSIDX, OBSVAL, LAMBDA) returns the N x 1 series
%   that best fits the values OBSVAL observed at the slice indices OBSIDX while
%   staying smooth. Smoothness is enforced with a second-difference penalty of
%   weight LAMBDA, so the result is a natural smoothing spline evaluated on the
%   slice grid: large LAMBDA gives a nearly straight line through the
%   observations, small LAMBDA follows them closely.
%
%   Unobserved slices are filled by the penalty alone, which is what lets a
%   handful of clicks define a centre line for the whole cord. The first and
%   last rows of the difference operator use a first difference instead, so the
%   series extrapolates linearly rather than curling at the ends.
%
%   YFIT = SMOOTHCORDSERIES(..., OBSWT) weights the observations, minimising
%   sum(OBSWT .* (YFIT(OBSIDX) - OBSVAL).^2) plus the penalty. Use it when some
%   observations are known to be worse than others: an observation of weight 0
%   is ignored entirely and its slice is filled in from its neighbours, which is
%   how AUTOCORDCENTERLINE stops sections with no measurable orientation from
%   dragging the fit around. Weights default to 1.
%
%   INPUTS
%     N       - number of slices (length of the output series).
%     obsidx  - indices of the observed slices.
%     obsval  - values observed at those slices, same length as OBSIDX.
%     lambda  - smoothness weight (>= 0).
%     obswt   - optional per-observation weight (>= 0), scalar or same length
%               as OBSIDX.
%
%   OUTPUT
%     yfit    - N x 1 double, the smoothed series.
%
%   See also AUTOCORDCENTERLINE, SPINAL_CORD_ALIGNER, FITCORDSERIES.

obsidx = obsidx(:);
obsval = obsval(:);

if nargin < 5 || isempty(obswt)
    obswt = ones(size(obsidx));
else
    obswt = obswt(:);
    if isscalar(obswt)
        obswt = repmat(obswt, size(obsidx));
    end
end

e  = ones(N, 1);
D2 = spdiags([e -2*e e], -1:1, N, N);
% linear (not curved) extrapolation at both ends
D2(1, :) = 0; D2(1, 1:2)   = [1 -1];
D2(N, :) = 0; D2(N, N-1:N) = [-1 1];

L = D2' * D2;
% accumulate, so a repeated slice index adds its weight instead of replacing it
W = sparse(obsidx, obsidx, obswt, N, N);
b = accumarray(obsidx, obswt .* obsval, [N 1]);

yfit = full((W + lambda * L) \ sparse(b));

end
