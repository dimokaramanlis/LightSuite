function yfit = fitCordSeries(N, obsidx, obsval, lambda)
%FITCORDSERIES Smoothing spline through observations, degenerate cases included.
%
%   YFIT = FITCORDSERIES(N, OBSIDX, OBSVAL, LAMBDA) is SMOOTHCORDSERIES with the
%   two cases it cannot solve handled:
%
%     no observations  - the series is undefined and comes back all NaN;
%     one observation  - the series is that constant. With a second-difference
%                        penalty alone the system is singular here, and a
%                        constant is the sensible limit: it happens when a
%                        single click has suppressed every prediction around it
%                        on a short cord.
%
%   See also SMOOTHCORDSERIES, MIXCORDOBSERVATIONS, SPINAL_CORD_ALIGNER.

if isempty(obsidx)
    yfit = nan(N, 1);
elseif isscalar(obsidx)
    yfit = repmat(obsval(1), N, 1);
else
    yfit = smoothCordSeries(N, obsidx, obsval, lambda);
end

end
