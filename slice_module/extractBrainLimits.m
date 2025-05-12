function [xrange, yrange] = extractBrainLimits(currslice, Nbuff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

currslice = single(currslice);
backval   = quantile(currslice(currslice>0), 0.01, 'all');
[ny, nx]  = size(currslice);
currslice = (currslice - backval)./backval;
minperc   = 0.25;
maxperc   = 0.75;

signaly   = single(max(currslice(:, round(nx*0.2):round(nx*0.8)), [], 2));
signaly   = signaly(:);
signalx   = single(max(currslice(round(ny*0.2):round(ny*0.8), :), [], 1));
signalx   = signalx(:);


minxsig   = signalx(1:round(nx*minperc));
xmin      = find(flip(minxsig)<2, 1, 'first'); %findfirst(flip(minxsig)<2);
if xmin > 0
    xmin      = numel(minxsig) - xmin;
else
    xmin = 1;
end

minysig   = signaly(1:round(ny*minperc));
ymin      = find(flip(minysig)<2,1,'first'); %findfirst(flip(minysig)<2);
if ymin > 0
    ymin      = numel(minysig) - ymin;
else
    ymin = 1;
end

maxxsig   = signalx(round(maxperc*nx):end);
xmax      = find(maxxsig<2, 1, 'first'); %findfirst(maxxsig<2);
if xmax > 0
    xmax      = round(nx*maxperc) + xmax;
else
    xmax = nx;
end

maxysig   = signaly(round(maxperc*ny):end);
ymax      = find(maxysig<2, 1, 'first'); %findfirst(maxysig<2);
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


yrange = yrange(1)-Nbuff:yrange(2)+Nbuff;
xrange = xrange(1)-Nbuff:xrange(2)+Nbuff;
yrange(yrange<1 | yrange > ny) = [];
xrange(xrange<1 | xrange > nx) = [];

end