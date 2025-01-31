function cf = plotAnnotationComparison(volume, anvol, dimplot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Ny    = size(volume, dimplot);
Nshow = 10;
ishow = round(linspace(0.15*Ny, 0.85*Ny, Nshow));

cf = figure;
cf.Position = [50 50 1700 900];
cf.Visible = 'off';
ppanel = panel();
ppanel.pack(2, Nshow/2)
ppanel.de.margin = 1;
ppanel.de.margintop = 6;
ppanel.margin = [1 1 1 8];


for ii = 1:Nshow
    islice = ishow(ii);
    [irow, icol] = ind2sub([2 Nshow/2], ii);
    ppanel(irow,icol).select();
    histim  = volumeIdtoImage(volume, [ islice dimplot]);
    atlasim = volumeIdtoImage(anvol, [ islice dimplot]);
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
    image(histim); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = gray;
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title(ishow(ii))
end


end