function cf = plotRegistrationComparison(slice, annotations, annstr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nannot = size(annotations, 3);
if ~isempty(annstr)
    assert(numel(annstr) == Nannot)
end


cf = figure;
cf.Position = [50 50 1500 600];
cf.Visible = 'off';
ppanel = panel();
ppanel.pack('h', Nannot)
ppanel.de.margin = 1;
ppanel.de.margintop = 6;
ppanel.margin = [1 1 1 8];
mksize = 2;

for ii = 1:Nannot
    ppanel(ii).select();
    atlasim = annotations(:, :, ii);
    atlasim = single(atlasim);
    
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
  
    image(slice); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = gray;
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',mksize)
    title(annstr{ii})
end


end