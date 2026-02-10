function [cinfo, cim] = cellDetectorNew(volumeuse, avgcellradius, ...
    sigmause, anisotropyratio, sdthres, bufferzone, saveimages)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
[Ny, Nx, Nz] = size(volumeuse);
usegpu       = isgpuarray(volumeuse);
%--------------------------------------------------------------------------
cubsize = ceil(4*sigmause.*anisotropyratio) + 1;
seuse   = strel('cuboid', cubsize);
%==========================================================================
% perform 3d bandpass filtering
% sdthres = [10 8];
[dff2, ~]  = spatial_bandpass_3d(volumeuse, avgcellradius, 3, 3, usegpu, anisotropyratio);
nhood      = floor(anisotropyratio*avgcellradius*10/2)*2+1;
backsd     = stdfilt3(dff2, nhood);
% background = volumeuse-dff2;
% background(background<0)=1e-3;
% %==========================================================================
% % calculate robust standard deviation to define noise level
% iy    = unique(round(linspace(1, Ny, round(20/anisotropyratio(1)))));
% ix    = unique(round(linspace(1, Nx, round(20/anisotropyratio(2)))));
% iz    = unique(round(linspace(1, Nz, round(20/anisotropyratio(3)))));
% 
% [xx, yy, zz] = meshgrid(ix, iy, iz);
% isamp = sub2ind([Ny, Nx, Nz], yy(:), xx(:), zz(:));
% currcheck = volumeuse(isamp);
% isamp(currcheck < 1) = [];
% 
% 
% warning off; 
% [idx, C ]  = kmeans(volumeuse(isamp), 2, 'Replicates', 10, 'MaxIter', 100);
% warning on; 
% [~, imax] = max(C);
% Npoints   = accumarray(idx, 1, [2 1], @sum);
% idsuse    = isamp;
% if min(Npoints)/max(Npoints)> 0.05
%     if max(C)/min(C) > 2
%         idsuse = isamp(idx==imax);
%     end
% end
% 
% dff2      = dff2 - median(dff2, 'all');
% backsd   = robustStd(dff2(idsuse));
% thresuse = median(volumeuse(idsuse));

% backsd    = robustStd(dff2(:))
%==========================================================================
ampsignal = dff2./backsd;
%==========================================================================
% detect local maxima and get candidate voxels
smax      = my_max(dff2, 4*sigmause, 1:3);
imgidx    = (dff2==smax) & (ampsignal > sdthres(1)); 
imgidx    = gather(imgidx);
imgidx    = imdilate(imgidx, seuse) & gather(ampsignal > sdthres(2));

cc    = bwconncomp(imgidx, 18);
cinfo = regionprops3(cc, gather(volumeuse),...
    'PrincipalAxisLength', 'EquivDiameter', 'VoxelList', ...
    'WeightedCentroid','MeanIntensity');
%==========================================================================
cim = [];
if size(cinfo, 1) > 0
    ccents = cinfo.WeightedCentroid;
    ikeepz = ccents(:, 3)> bufferzone(3,1) & (ccents(:, 3) < (Nz - bufferzone(3,2)));
    ikeepx = ccents(:, 1)> bufferzone(1,1) & (ccents(:, 1) < (Nx - bufferzone(1,2)));
    ikeepy = ccents(:, 2)> bufferzone(2,1) & (ccents(:, 2) < (Ny - bufferzone(2,2)));
    
    elips    = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
    celldiam = cinfo.EquivDiameter * prod(1./anisotropyratio)^(1/3);
    ilong    = elips>2.5;
    ismall   = celldiam<avgcellradius;

    % goodintensity = cinfo.MeanIntensity(~ilong & ~ismall);
    % if numel(goodintensity) > 100
    %     intthres = quantile(goodintensity,0.01);
    % else
    %     intthres = 0;
    % end
    intthres = 0;
    ihigh = cinfo.MeanIntensity > intthres;
    % %%
    % [~, bb]= pca(zscore([elips,celldiam,cinfo.MeanIntensity]));
    % clf; hold on;
    % [idxcells, C] = kmedoids(bb(:,1:2), 4, 'Replicates',5);
    % for ii =  1:4
    %     plot(bb(idxcells==ii,1),bb(idxcells==ii,2),'o')
    % end
    %%
    ikeep  = ikeepz & ikeepx & ikeepy & ~ilong & ~ismall & ihigh;
    cinfo  = cinfo(ikeep, :);
    if saveimages
        cim = getCellImages2D(dff2, cinfo, sigmause*6);
    end
end
%=========================================================================
% %intermediate plotting
% if size(cinfo, 1) > 0
%     allvoxels =  cat(1,cinfo.VoxelList{:});
%     indtest   = sub2ind(size(dff2),allvoxels(:,2), allvoxels(:,1), allvoxels(:,3));
%     imgidx = false(size(imgidx));
%     imgidx(indtest) = true;
%     imshowpair(uint8(max(ampsignal,[],3)*255/sdthres(1)),max(imgidx,[],3));
% end
%=========================================================================   
end

