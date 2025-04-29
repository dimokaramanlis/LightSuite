function [cpoints, prem] = cellDetector(volumeuse, avgcellradius, sigmause, anisotropyratio, thresSNR)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz] = size(volumeuse);
%--------------------------------------------------------------------------
cubsize = 4*sigmause ;
seuse   = strel('cuboid', cubsize);

usegpu  = isgpuarray(volumeuse);

% % let's learn cells from data
% 4 times as big as the cell
% 1.5 times smaller than the cell
[fout, ~] = spatial_bandpass_3d(volumeuse, avgcellradius, 4, 1.5, usegpu, anisotropyratio);

% find maxima in filtered space
smin           = -my_min(-fout, 2*sigmause, 1:3);

% we keep maxima with a minimum snr
% this snr should be probably in filtered space
imgidx = (fout>smin-1e-3) & (fout > thresSNR(1));
imgidx = gather(imgidx);
imgidx = imdilate(imgidx, seuse) & gather(fout > thresSNR(2));

% which properties are needed???
cinfo  = regionprops3(imgidx, gather(volumeuse),...
    'PrincipalAxisLength', 'EquivDiameter', 'VoxelList', 'WeightedCentroid','MeanIntensity');
%--------------------------------------------------------------------------
% for debugging

% for ii = 1:size(volumeuse,3)
%     imshowpair(uint8(volumeuse(:,:,ii)*255/thresSNR(1)),imgidx(:,:,ii));
%     pause; 
% end
%--------------------------------------------------------------------------
% improve equivalent diameter

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
    elips   = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
    ilong   = elips>2.5;
    ismall  = cinfo.EquivDiameter<avgcellradius/2;
    ilow    = cinfo.MeanIntensity < thresSNR(2);
    
    iweird = find(ilong | ismall | ilow) ;
    prem   = numel(iweird)/size(cinfo, 1);
    
    % allvoxels =  cat(1,cinfo(iweird,:).VoxelList{:});
    % indtest = sub2ind(size(fout),allvoxels(:,1), allvoxels(:,2), allvoxels(:,3));
    % imgidx = false(size(imgidx));
    % imgidx(indtest) = true;
    % imshowpair(uint8(max(volumeuse,[],3)*255/thresSNR(1)),max(imgidx,[],3));

    cinfo(iweird, :) = [];    
end
%--------------------------------------------------------------------------
% package cell properties
cpoints = [cinfo.WeightedCentroid cinfo.EquivDiameter cinfo.MeanIntensity];
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