function ptcloud = extractHighSFVolumePoints(voluse, pxsize, varargin)
%EXTRACTMATCHINGGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz] = size(voluse);

if nargin < 3
    rykeep = [1 Ny];
else
    rykeep = varargin{1};
end

scalefilter = 100/pxsize;
imhigh     = spatial_bandpass_3d(voluse, scalefilter, 3, 3, true);
ptthres    = quantile(imhigh, 0.99, 'all')/2;
ipts       = find(imhigh>ptthres);
[rr,cc,dd] = ind2sub(size(voluse), ipts);
X          = [cc,rr,dd];
ikeep      = rr > rykeep(1) & rr<rykeep(2);
X          = X(ikeep, :);
ptcloud    = pointCloud(X);
Npts       = nnz(ikeep);
ptcloud    = pcdownsample(ptcloud,'random', 10000/Npts, 'PreserveStructure',true);

 % scatter3(ptcloud.Location(:,1),ptcloud.Location(:,2),ptcloud.Location(:,3),2)
end

