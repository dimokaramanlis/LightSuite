function cf = plotAnnotationComparison(volume, anvol, dimplot, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%   plotAnnotationComparison(volume, anvol, dimplot, pxsize, ...)
%   Optional name-value pairs:
%       'MarkerSize' - size of the boundary markers (default 1)
%       'Visible'    - show the figure (true/'on') or keep it hidden
%                      (false/'off', default)

p = inputParser;
p.addOptional('pxsize', [1 1 1], @(x) isnumeric(x) && isvector(x));
p.addParameter('MarkerSize', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Visible', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)) || ...
    any(strcmpi(x, {'on','off'})));
p.parse(varargin{:});

pxsize = p.Results.pxsize;
mksize = p.Results.MarkerSize;
if ischar(p.Results.Visible) || isstring(p.Results.Visible)
    visible = char(p.Results.Visible);
else
    visstates = {'off','on'};
    visible = visstates{logical(p.Results.Visible) + 1};
end

pxsize(dimplot) = [];

Ny    = size(volume, dimplot);
Nshow = 8;
ishow = round(linspace(0.15*Ny, 0.85*Ny, Nshow));

cf = figure;
cf.Position = [50 50 1700 900];
cf.Visible = visible;
ppanel = panel();
ppanel.pack(2, Nshow/2)
ppanel.de.margin = 1;
ppanel.de.margintop = 6;
ppanel.margin = [1 1 1 8];

atlassize = size(anvol);
atlassize(dimplot) = [];
xvals     = 1:pxsize(2):atlassize(2)*pxsize(2);
yvals     = 1:pxsize(1):atlassize(1)*pxsize(1);

for ii = 1:Nshow
    islice = ishow(ii);
    [irow, icol] = ind2sub([2 Nshow/2], ii);
    ppanel(irow,icol).select();
    histim  = volumeIdtoImage(volume, [ islice dimplot]);
    atlasim = volumeIdtoImage(anvol, [ islice dimplot]);
    atlasim = single(atlasim);
    
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
  
    image(xvals, yvals, histim); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = gray;
    line(col*pxsize(2), row*pxsize(1), 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',mksize)
    title(ishow(ii))
end


end