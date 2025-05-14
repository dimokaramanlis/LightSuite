function [cpoints, prem] = cellDetector2D(sliceuse, avgcellradius, sigmause, thresSNR)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx] = size(sliceuse);
%--------------------------------------------------------------------------
cubsize = 4*sigmause ;
seuse   = strel('rectangle', cubsize);

usegpu  = isgpuarray(sliceuse);

% % let's learn cells from data
% 4 times as big as the cell
% 1.5 times smaller than the cell
[fout, ~] = spatial_bandpass(sliceuse, avgcellradius, 4, 1.5, usegpu);

% find maxima in filtered space
smin           = -my_min(-fout, 2*sigmause, 1:2);

% we keep maxima with a minimum snr
% this snr should be probably in filtered space
imgidx = (fout>smin-1e-3) & (fout > thresSNR(1));
imgidx = gather(imgidx);
imgidx = imdilate(imgidx, seuse) & gather(fout > thresSNR(2));
% imshowpair(uint8(sliceuse*255/thresSNR(1)),imgidx);
% which properties are needed???
cinfo  = regionprops(imgidx, gather(sliceuse),...
    'Circularity', 'EquivDiameter', 'PixelList', 'WeightedCentroid','MeanIntensity');
%--------------------------------------------------------------------------
% for debugging

%-------------------------------------------------------------------------
% improve equivalent diameter and cell filtering



% dty = -sigmause(1)*4:sigmause(1)*4;
% dtx = -sigmause(2)*4:sigmause(2)*4;
% 
% X = nan(size(cinfo,1), numel(dtx)*numel(dty),'single');
% for icell = 1:size(cinfo,1)
%     ccent  =  round(cinfo(icell).WeightedCentroid);
%     xind = ccent(1) + dtx;
%     yind = ccent(2) + dty;
%     if any(xind<1)|any(yind<1)|any(xind>Nx)|any(yind>Ny)
%         continue
%     end
%     currvol = sliceuse(yind, xind);
%     X(icell,:) = reshape(currvol, 1, []);
% end



% [cellimages, cellproperties] = getCellImages2D(volumeuse, cinfo)

% we find unique pairs of x-y pixels, take their sum as the area of the
% cell, then find diameter. to remove z-blurring
%--------------------------------------------------------------------------

% remove weird cells: elongated or small or with low intensity

% CAN WE DO BETTER?


% alldist = squareform(pdist(cinfo.Centroid));
% alldists = sort(alldist, 2, 'ascend');
% localdist = median(alldists(:, 2:6), 2);
% 
% [localfun] = histc(alldists,linspace(0,400,20),2);
% localfun = localfun(:,1:end-1)./sum(localfun(:,1:end-1),2);
% [aa,bb] = pca(localfun);
% iweird = find(localdist>80);
prem   = 0;
if ~isempty(cinfo)
    ilong   = [cinfo(:).Circularity]<0.7;
    ismall  = [cinfo(:).EquivDiameter]<avgcellradius/2;
    ilow    = [cinfo(:).MeanIntensity] < thresSNR(2);
    
    iweird = find(ilong | ismall | ilow) ;
    prem   = numel(iweird)/size(cinfo, 1);
    
    % allvoxels =  cat(1,cinfo(iweird).PixelList);
    % indtest = sub2ind(size(fout),allvoxels(:,2), allvoxels(:,1));
    % imgidx = false(size(imgidx));
    % imgidx(indtest) = true;
    % imshowpair(uint8(sliceuse*255/thresSNR(1)),imgidx);

    cinfo(iweird) = [];    
end
%--------------------------------------------------------------------------
% package cell properties
cpoints = [cat(1,cinfo.WeightedCentroid) cat(1,cinfo.EquivDiameter) cat(1,cinfo.MeanIntensity)];
%--------------------------------------------------------------------------
   
end



% % 
% 
% dty = -sigmause(1)*4:sigmause(1)*4;
% dtx = -sigmause(2)*4:sigmause(2)*4;
% dtz = -sigmause(3)*4:sigmause(3)*4;
% 
% X = nan(size(cinfo,1), numel(dtx)*numel(dty)*numel(dtz),'single');
% for icell = 1:size(cinfo,1)
%     ccent  =  round(cinfo.WeightedCentroid(icell,:));
%     xind = ccent(1) + dtx;
%     yind = ccent(2) + dty;
%     zind = ccent(3) + dtz;
%     if any(zind<1)|any(xind<1)|any(yind<1)|any(xind>Nx)|any(yind>Ny)|any(zind>Nz)
%         continue
%     end
%     currvol = volumeuse(yind, xind, zind);
%     X(icell,:) = reshape(currvol, 1, []);
% end
% 
% irem = any(isnan(X),2);
% XX = X(~irem, :);
% XX = XX./sqrt(sum(XX.^2,2));
% [aa,bb] = pca(XX);
% aplot = reshape(aa(:,2),numel(dty),numel(dtx),numel(dtz));
% explot = reshape(XX(isort(10),:),numel(dty),numel(dtx),numel(dtz));
% 
% idxgroup = aa(:,2)>0;