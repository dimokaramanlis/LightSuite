function cf = visualizeTransformMatch(commonvol, matchvols, dimplot, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

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

Ntypes  = size(matchvols, 4);
Ny    = size(commonvol, dimplot);
Nshow = 5;
ishow = round(linspace(0.15*Ny, 0.85*Ny, Nshow));

cf = figure;
cf.Position = [50 50 1750 250*Ntypes];
cf.Visible = visible;
ppanel = panel();
ppanel.pack(Ntypes, Nshow)
ppanel.de.margin = 1;
ppanel.de.margintop = 4;
ppanel.margin = [1 1 1 8];

atlassize = size(commonvol);
atlassize(dimplot) = [];
xvals     = 1:pxsize(2):atlassize(2)*pxsize(2);
yvals     = 1:pxsize(1):atlassize(1)*pxsize(1);

for ii = 1:Nshow
    islice = ishow(ii);
    commonimg = volumeIdtoImage(commonvol, [ islice dimplot]);

    for itype = 1:Ntypes
        histim  = volumeIdtoImage(matchvols(:, :, :, itype), [ islice dimplot]);
        C = imfuse(commonimg,histim,'falsecolor','Scaling','joint');
        axcurr = ppanel(itype, ii).select();
        image(C);
        axis equal; axis tight;
        axcurr.YDir = 'reverse'; axcurr.Visible = 'off';
        axcurr.Title.Visible = 'on';
        if itype == 1
            title(ishow(ii))
        end
    end
    % histim  = volumeIdtoImage(volume, [ islice dimplot]);
    % atlasim = volumeIdtoImage(anvol, [ islice dimplot]);
    % atlasim = single(atlasim);
    % 
    % av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    % [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
    % 
    % image(xvals, yvals, histim); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    % ax.Title.Visible = 'on';
    % line(col*pxsize(2), row*pxsize(1), 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',mksize)
    % title(ishow(ii))
end



end