function [xrange, yrange] = extractBrainLimits3(currslice, Nbuff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Nwidth = floor(Nbuff/4)*2 + 1;
[ny, nx]  = size(currslice);


% matuse   = ones(Nwidth);
% Nmed     = floor(sum(matuse, 'all')*0.5);
% filtim = ordfilt2(currslice, Nmed, matuse);
% dfimg  = (single(currslice)- single(filtim))./single(filtim);
% dfimg(isinf(dfimg)) = 0;

backval = single(quantile(currslice(currslice>0), 0.01, 'all'));
dfimg  = (single(currslice)- backval)./backval;
dfimg(isinf(dfimg) | dfimg < 0) = 0;


minperc   = 0.25;
maxperc   = 0.75;
xlook = round(nx*minperc):round(nx*maxperc);
ylook = round(ny*minperc):round(ny*maxperc);
centvals = dfimg(ylook, xlook);
thresuse = max(0.5, quantile(centvals,0.99,'all')/4);
ipx = find(dfimg(:) > thresuse);
[row, col] = ind2sub(size(dfimg), ipx);
iremy = row < Nbuff | row > ny - Nbuff;
iremx = col < Nbuff | col > nx - Nbuff;
irem = iremy | iremx;

pccloud = pointCloud([col(~irem), row(~irem), zeros(size(row(~irem)))]);
Npts    = pccloud.Count;
targetnum = min(2e4, Npts);
pccloud = pcdownsample(pccloud, 'random', targetnum/Npts);
pccloud = pcdenoise(pccloud, 'NumNeighbors',100);

% gmm = fitgmdist(pccloud.Location(:, 1:2),1);
aa = pcsegdist(pccloud, Nbuff/2);
clustersizes = accumarray(aa, 1, [], @sum);
clustkeep = find(clustersizes/sum(clustersizes) > 0.1);
clustids = ismembc(aa, uint32(clustkeep));

% clf;
% imagesc(dfimg, [0 thresuse]); hold on;
% plot(pccloud.Location(clustids, 1),pccloud.Location(clustids, 2),'r.')
%%
xmin = min(pccloud.Location(clustids, 1));
xmax = max(pccloud.Location(clustids, 1));
ymin = min(pccloud.Location(clustids, 2));
ymax = max(pccloud.Location(clustids, 2));

xrange  = [xmin xmax];
yrange  = [ymin ymax];
% if (xmax-xmin)/nx < 0.5
%     Nbuffx = 2*Nbuff;
% else
%     Nbuffx = Nbuff;
% end
% if (ymax-ymin)/ny < 0.5
%     Nbuffy = 2*Nbuff;
% else
%     Nbuffy = Nbuff;
% end
Nbuffy = Nbuff;
Nbuffx = Nbuff;

yrange = yrange(1)-Nbuffy:yrange(2)+Nbuffy;
xrange = xrange(1)-Nbuffx:xrange(2)+Nbuffx;
yrange(yrange<1 | yrange > ny) = [];
xrange(xrange<1 | xrange > nx) = [];

end