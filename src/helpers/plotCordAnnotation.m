function cf = plotCordAnnotation(volume, anvol)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



atlassize = size(anvol);
Nshow = 12;
Nlong = 4;
anncol = [1 0.8 0.5];
cf = figure;
cf.Position = [50 50 1500 850];
cf.Visible = 'off';
%%
ppanel = panel();
ppanel.pack('h', {0.5 0.5});
ppanel(1).pack('h', 2)
for ii = 1:2
    ppanel(1, ii).pack('h', Nlong);
end
ppanel(2).pack(Nshow/3, 3)

ppanel.de.margin = 1;
ppanel(1,2).marginleft = 8;
ppanel(2).marginleft = 8;
% ppanel(1).de.margintop = 6;
ppanel.margin = [0.5 1 0.5 6];
%--------------------------------------------------------------------------
Nlength    = atlassize(3);
ishowslice = round(linspace(0.02*Nlength, 0.98*Nlength, Nshow));

for ii = 1:Nshow
    islice = ishowslice(ii);
    [irow, icol] = ind2sub([Nshow/3 3], ii);
    ppanel(2, irow,icol).select();
    histim  = volume(:,:,islice);
    atlasim = anvol(:,:,islice);
    atlasim = single(atlasim);
    
    % av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    % [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));

     av_warp_boundaries = round(conv2(atlasim,ones(3)./9,'same')) ~= atlasim;
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
  
    image(histim); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    ax.YDir = 'reverse'; ax.Colormap = gray;
    mksize = 1;
    line(col, row, 'Marker','.','LineStyle','none', 'Color',anncol,'MarkerSize',mksize)
    currsize= size(atlasim);
    text(currsize(2)*0.95, currsize(1)*0.95, sprintf('%d', ishowslice(ii)), 'Color', 'w',...
        'HorizontalAlignment','right')
end
%--------------------------------------------------------------------------

Nlength1   = atlassize(1);
ishowslice = round(linspace(0.25*Nlength1, 0.75*Nlength1, Nlong));

for ii = 1:Nlong
    islice = ishowslice(ii);
    ppanel(1,1,ii).select(); hold on;
    histim  = squeeze(volume(islice,:,:))';
    atlasim = squeeze(anvol(islice,:,:))';
    % atlasim = single(atlasim);
    % atlasim = squeeze(avbound(islice,:,:))';
    % av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    % [row,col] = ind2sub(size(atlasim), find(atlasim));

    av_warp_boundaries = round(conv2(atlasim,ones(3)./9,'same')) ~= atlasim;
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));

  
    image(histim); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    ax.YDir = 'reverse'; ax.Colormap = gray; ax.Title.Visible = 'on';
    % mksize = 1;
    % line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',mksize)
    scatter(col, row, 1, 'filled', 'MarkerFaceColor', anncol)  
    title(ishowslice(ii))
end
%--------------------------------------------------------------------------


Nlength2   = atlassize(2);
ishowslice = round(linspace(0.25*Nlength2, 0.75*Nlength2, Nlong));

for ii = 1:Nlong
    islice = ishowslice(ii);
    ppanel(1,2,ii).select(); hold on;
    histim  = squeeze(volume(:,islice,:))';
    atlasim = squeeze(anvol(:,islice,:))';
    av_warp_boundaries = round(conv2(atlasim,ones(3)./9,'same')) ~= atlasim;
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    image(histim); ax = gca; ax.Visible = 'off';axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = gray;
    mksize = 0.2;
    % line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',mksize)
    scatter(col, row, 1, 'filled', 'MarkerFaceColor', anncol)
    title(ishowslice(ii))
end
%--------------------------------------------------------------------------
%%
end