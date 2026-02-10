function mousemap = atlasCellMap(mousecells, avsize, sigma, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if iscell(mousecells)
    Nmice     = numel(mousecells);
    Ntot      = sum(cellfun(@(x) size(x,1), mousecells));
    % pointlocs = cat(1, mousecells{:});
else
    Nmice = 1;
    mousecells = {mousecells};
    Ntot  = size(pointlocs, 1);
end

if nargin < 4
    wts = ones(Nmice, 1)/Nmice;
else
    wts = varargin{1}(:);
    wts = wts/sum(wts);
end

avplot       = ndSparse.spalloc(avsize, Ntot);
for ii = 1:Nmice
    avplot       = avplot + ndSparse.build(double(mousecells{ii}(:,[2 1 3])), wts(ii),avsize);
end

% pointlocs    = double(pointlocs);
avplot       = single(full(avplot));

Nmid         = avsize(3)/2;
avplot       = (avplot(:,:, 1:Nmid) + flip(avplot(:,:, (Nmid + 1):end), 3))/2;
mousemap     = imgaussfilt3(avplot, sigma);

end