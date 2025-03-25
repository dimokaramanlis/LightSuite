function mousemap = atlasCellMap(mousecells,avsize, sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if iscell(mousecells)
    Nmice     = numel(mousecells);
    pointlocs = cat(1, mousecells{:});
else
    Nmice = 1;
    pointlocs = mousecells;
end

avplot       = ndSparse.build(pointlocs(:,[2 1 3]), 1,avsize);
avplot       = single(full(avplot));

Nmid         =  avsize(3)/2;
avplot       = (avplot(:,:, 1:Nmid) + flip(avplot(:,:, (Nmid + 1):end), 3))/2;
mousemap     = imgaussfilt3(avplot, sigma)/Nmice;

end