function [xrange, yrange] = extractBrainLimits2(currslice, Nbuff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Nwidth = floor(Nbuff/4)*2 + 1;
[ny, nx]  = size(currslice);

% filtim = medfilt2(currslice, [Nwidth Nwidth]);

matuse   = ones(Nwidth);
Nmed     = floor(sum(matuse, 'all')*0.5);
filtim = ordfilt2(currslice, Nmed, matuse);

dfimg  = (single(currslice)- single(filtim))./single(filtim);
dfimg(isinf(dfimg)) = 0;
%%
ipx = find(dfimg(:) > 0.2);
[row, col] = ind2sub(size(dfimg), ipx);
iremy = row < Nbuff | row > ny - Nbuff;
iremx = col < Nbuff | col > nx - Nbuff;
irem = iremy | iremx;

pccloud = pointCloud([col(~irem), row(~irem), zeros(size(row(~irem)))]);
pccloud = pcdenoise(pccloud, 'NumNeighbors',100);
Npts    = pccloud.Count;
pccloud = pcdownsample(pccloud, 'random', 1e4/Npts);
% gmm = fitgmdist(pccloud.Location(:, 1:2),1);
aa = pcsegdist(pccloud, 100);
clustersizes = accumarray(aa, 1, [], @sum);
clustkeep = find(clustersizes/sum(clustersizes) > 0.3);
clustids = ismembc(aa, uint32(clustkeep));
% gauss.mu = gmm.mu';
% gauss.sigma = gmm.Sigma;
% c = getEllipse(gauss, 3);
clf;
imagesc(dfimg, [0 1]); hold on;
plot(pccloud.Location(clustids, 1),pccloud.Location(clustids, 2),'r.')
%%
outim = single(currslice);
outim = outim - (outim*single(xbase)')*single(xbase)/norm(single(xbase))^2;
outim = outim - single(ybase)*(outim'*single(ybase))'/norm(single(ybase))^2;


%%
% dfimg  = dfimg.*(dfimg>0);
% dfimg  = dfimg/quantile(dfimg, 0.99,'all');
% 
% sliceinfo.px_process = 1.25;
% scalefilter = 100/sliceinfo.px_process;
% singslice = single(currslice);
% imhigh   = spatial_bandpass(singslice, scalefilter, Inf, 1, true);

% backval   = quantile(currslice(currslice>0), 0.01, 'all');
[ny, nx]  = size(currslice);
% currslice = (currslice - backval)./backval;
minperc   = 0.25;
maxperc   = 0.75;
signalt   = 0.1; % background threshold

signaly   = single(quantile(dfimg(:, round(nx*minperc):round(nx*maxperc)), 0.99, 2));
signaly   = movmedian(signaly(:), round(Nbuff/4));
signalx   = single(quantile(dfimg(round(ny*minperc):round(ny*maxperc), :), 0.99, 1));
signalx   = movmedian(signalx(:), round(Nbuff/4));


% signaly   = single(max(currslice(:, round(nx*0.2):round(nx*0.8)), [], 2));
% signaly   = movmedian(signaly(:), round(Nbuff/4));
% signalx   = single(max(currslice(round(ny*0.2):round(ny*0.8), :), [], 1));
% signalx   = movmedian(signalx(:), round(Nbuff/4));


minxsig   = signalx(1:round(nx*minperc));
xmin      = find(flip(minxsig)<signalt, 1, 'first'); %findfirst(flip(minxsig)<2);
if xmin > 0
    xmin      = numel(minxsig) - xmin;
else
    xmin = 1;
end

minysig   = signaly(1:round(ny*minperc));
ymin      = find(flip(minysig)<signalt,1,'first'); %findfirst(flip(minysig)<2);
if ymin > 0
    ymin      = numel(minysig) - ymin;
else
    ymin = 1;
end

maxxsig   = signalx(round(maxperc*nx):end);
xmax      = find(maxxsig<signalt, 1, 'first'); %findfirst(maxxsig<2);
if xmax > 0
    xmax      = round(nx*maxperc) + xmax;
else
    xmax = nx;
end

maxysig   = signaly(round(maxperc*ny):end);
ymax      = find(maxysig<signalt, 1, 'first'); %findfirst(maxysig<2);
if ymax > 0
    ymax      = round(ny*maxperc) + ymax;
else
    ymax = ny;
end


% 
% 
% xmin = findchangepts(signalx(1:round(nx*0.2)), 'MaxNumChanges',1);
% if isempty(xmin)
%     xmin = findchangepts(signalx(1:round(nx*0.2)), 'MaxNumChanges',2);
%     if ~isempty(xmin)
%         xmin = xmin(2);
%     else
%         xmin = 1;
%     end
% end
% 
% xmax      = findchangepts(signalx(round(nx*0.8):end), 'MaxNumChanges',1);
% if isempty(xmax) 
%     xmax = findchangepts(signalx(round(nx*0.8):end), 'MaxNumChanges',2);
%     if ~isempty(xmax)
%         xmax = xmax(1);
%     else
%         xmax = nx;
%     end
% end
% 
% ymin = findchangepts(signaly(1:round(ny*0.2)), 'MaxNumChanges',1);
% if isempty(ymin)
%     ymin = findchangepts(signaly(1:round(ny*0.2)), 'MaxNumChanges',2);
%     if ~isempty(ymin)
%         ymin = ymin(2);
%     else
%         ymin = 1;
%     end
% end
% 
% ymax = findchangepts(signaly(round(ny*0.8):end), 'MaxNumChanges',1);
% if isempty(ymax)
%     ymax = findchangepts(signaly(round(ny*0.8):end), 'MaxNumChanges',2);
%     if ~isempty(ymax)
%         ymax = ymax(1);
%     else
%         ymax = ny;
%     end
% end
% 

xrange  = [xmin xmax];
yrange  = [ymin ymax];
if (xmax-xmin)/nx < 0.6
    Nbuffx = 2*Nbuff;
else
    Nbuffx = Nbuff;
end
if (ymax-ymin)/ny < 0.6
    Nbuffy = 2*Nbuff;
else
    Nbuffy = Nbuff;
end

yrange = yrange(1)-Nbuffy:yrange(2)+Nbuffy;
xrange = xrange(1)-Nbuffx:xrange(2)+Nbuffx;
yrange(yrange<1 | yrange > ny) = [];
xrange(xrange<1 | xrange > nx) = [];

end