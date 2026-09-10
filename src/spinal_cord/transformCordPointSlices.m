function ptsuse = transformCordPointSlices(ptsuse, tforms)
%TRANSFORMCORDPOINTSLICES Apply the per-slice cord transforms to points.
%
%   PTSOUT = TRANSFORMCORDPOINTSLICES(PTS, TFORMS) transforms the [x y] columns
%   of PTS with the transform of the slice its rounded z falls on. It is the
%   nearest-slice counterpart of TRANSFORMCORDIMAGESLICES and is meant for point
%   clouds sampled on slice centres.
%
%   Detected cells do not sit on slice centres, so CORDPOINTSTOATLAS uses
%   CORDSTRAIGHTENPOINTS instead, which interpolates the centre and the angle at
%   the point's own fractional z.
%
%   See also TRANSFORMCORDIMAGESLICES, CORDSTRAIGHTENPOINTS.

Nslices   = numel(tforms);
sliceinds = accumarray(ptsuse(:, 3), (1:size(ptsuse, 1))', [Nslices 1], @(x) {x});

for islice = 1:Nslices
    currinds = sliceinds{islice};
    if numel(currinds) > 0
        currpts = ptsuse(currinds, 1:2);
        ptsuse(currinds, 1:2) = tforms(islice).transformPointsForward(currpts);
    end
end

end
