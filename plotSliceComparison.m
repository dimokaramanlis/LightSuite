function cf = plotSliceComparison(volume,atlas, dimplot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Ny    = size(volume, dimplot);
Nshow = 10;
ishow = round(linspace(0.15*Ny, 0.85*Ny, Nshow));

cf = figure;
ppanel = panel();
ppanel.pack(2, Nshow/2)
ppanel.de.margin = 1;
% ppanel.margin = [1 1 1 1];

for ii = 1:Nshow
    islice = ishow(ii);
    [irow, icol] = ind2sub([2 Nshow/2], ii);
    ppanel(irow,icol).select();
    histim  = volumeIdtoImage(volume, [ islice dimplot]);
    atlasim = volumeIdtoImage(atlas, [ islice dimplot]);
    imshowpair(histim, atlasim)
    title(ishow(ii))
end


end