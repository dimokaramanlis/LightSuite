function ptvol = spinalCordPointCloud(cordvol, ikeep, volgood)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%-------------------------------------------------------------------------
% ymax     = max(cordvol,[],'all');
% T        = ymax * graythresh(cordvol);
% randvals = getRandomValsFromArray(cordvol(volgood), 2e4);
% tthres   = quantile(randvals, 0.5);
%-------------------------------------------------------------------------
Nslices  = size(cordvol, 3);
ptlist   = cell(Nslices, 1);
thresuse = 1;
for islice = 1:Nslices
    scurr = single(medfilt2(cordvol(:, :,islice)));
    igood = volgood(:, :, islice);

    % tthres = quantile(scurr(igood),0.5);
    % ipts = find(scurr > tthres);
    % % 
    gcurr = imgradient(scurr);
    gout  = gcurr./scurr;

    ipts  = find(gout>thresuse & igood);
    if numel(ipts) > 0
        [rr,cc]        = ind2sub(size(scurr), ipts);
        ptslice        = [cc, rr,  islice*ones(size(rr))];
        ptcurr         = pointCloud(ptslice);
        ptcurr         = pcdownsample(ptcurr, 'nonuniformGridSample', 12);
        ptlist{islice} = ptcurr.Location;    
    end
    % imagesc(scurr>quantile(scurr(igood),0.5))
    % pause
end
ptvol          = single(cat(1, ptlist{:}));
irem           = ptvol(:, 3) < ikeep(1) | ptvol(:, 3) > ikeep(2);
ptvol(irem, :) = [];

% 
% tvpts    = tvpts(:, [2 1 3]);
% 
% 
% ls_cloud       = extractVolumePointsGradient(single(cordvol), 1, 5);
% Nptsfin        = max(6, floor(ls_cloud.Count/1e5));
% ls_cloud_use   = pcdownsample(ls_cloud, 'nonuniformGridSample', Nptsfin);
% ls_cloud_use   = pcdenoise(ls_cloud_use);
% ptvol          = single(ls_cloud_use.Location);


end