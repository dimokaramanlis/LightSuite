function [dout, dtext] = compareDivisionDecoding(sin, ytopred, normfac, ikeep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

densitiesuse      = sin.counts./sin.volumes;
densitiesuse      = densitiesuse./normfac';
[divnames, ~, ic] = unique(sin.names(:,3));
decodedivisions   = nan(numel(divnames), 1);

for iname = 1:numel(divnames)
    currdata = densitiesuse(ic==iname, :);
    decodedivisions(iname) = decodeTaskFromSignal(currdata(:, ikeep), ytopred(ikeep));
end
dall = decodeTaskFromSignal(densitiesuse(:, ikeep), ytopred(ikeep));
dout = [dall; decodedivisions];
dtext = ['All'; divnames];
end