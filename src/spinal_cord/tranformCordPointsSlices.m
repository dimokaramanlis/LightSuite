function ptsuse = tranformCordPointsSlices(ptsuse, tforms)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


Nslices   = numel(tforms);
sliceinds = accumarray(ptsuse(:, 3), (1:size(ptsuse, 1))', [Nslices 1], @(x) {x});

for islice = 1:numel(tforms)
    currinds = sliceinds{islice};
    if numel(currinds) > 0
        currpts = ptsuse(currinds, 1:2);
        ptsuse(currinds, 1:2) = tforms(islice).transformPointsForward(currpts);
    end
end

end