function X = extractHighSFVolumePoints(voluse, pxsize, varargin)
%EXTRACTMATCHINGGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz] = size(voluse);
rykeep       = [1 Ny];

if nargin < 3
    usegpu = true;
else
    usegpu = varargin{1};
end

scalefilter = 100/pxsize;
imhigh     = spatial_bandpass_3d(voluse, scalefilter, 3, 3, usegpu);
ptthres    = quantile(imhigh, 0.99, 'all')/2;
ipts       = find(imhigh>ptthres);
[rr,cc,dd] = ind2sub(size(voluse), ipts);
X          = [cc,rr,dd];
ikeep      = rr > rykeep(1) & rr<rykeep(2);
X          = X(ikeep, :);
pcdown     = pcdownsample(pointCloud(X), 'random', 0.1);
X          = gather(pcdown.Location);
 % scatter3(ptcloud.Location(:,1),ptcloud.Location(:,2),ptcloud.Location(:,3),2)
end

