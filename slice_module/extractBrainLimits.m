function [xrange, yrange] = extractBrainLimits(currslice)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[ny, nx] = size(currslice);
signaly = single(max(currslice(:, round(nx*0.2):round(nx*0.8)), [], 2));
signalx = single(max(currslice(round(ny*0.2):round(ny*0.8), :), [], 1));

xmin = findchangepts(signalx(1:round(nx*0.2)), 'MaxNumChanges',1);
if isempty(xmin)
    xmin = findchangepts(signalx(1:round(nx*0.2)), 'MaxNumChanges',2);
    if ~isempty(xmin)
        xmin = xmin(2);
    else
        xmin = 1;
    end
end

xmax      = findchangepts(signalx(round(nx*0.8):end), 'MaxNumChanges',1);
if isempty(xmax) 
    xmax = findchangepts(signalx(round(nx*0.8):end), 'MaxNumChanges',2);
    if ~isempty(xmax)
        xmax = xmax(1);
    else
        xmax = nx;
    end
end

ymin = findchangepts(signaly(1:round(ny*0.2)), 'MaxNumChanges',1);
if isempty(ymin)
    ymin = findchangepts(signaly(1:round(ny*0.2)), 'MaxNumChanges',2);
    if ~isempty(ymin)
        ymin = ymin(2);
    else
        ymin = 1;
    end
end

ymax = findchangepts(signaly(round(ny*0.8):end), 'MaxNumChanges',1);
if isempty(ymax)
    ymax = findchangepts(signaly(round(ny*0.8):end), 'MaxNumChanges',2);
    if ~isempty(ymax)
        ymax = ymax(1);
    else
        ymax = ny;
    end
end


xrange  = [xmin round(nx*0.8) + xmax];
yrange  = [ymin round(ny*0.8)+ymax];
end