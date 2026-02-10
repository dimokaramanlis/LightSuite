function [cinfo, cim] = cellDetectorNew2(volumeuse, avgcellradius, ...
    sigmause, anisotropyratio, sdthres, bufferzone, saveimages)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
[Ny, Nx, Nz] = size(volumeuse);
usegpu       = isgpuarray(volumeuse);
sdthres      = [1 2];
%--------------------------------------------------------------------------
cubsize = ceil(4*sigmause.*anisotropyratio) + 1;
seuse   = strel('cuboid', cubsize);
%==========================================================================
% perform 3d bandpass filtering
% sdthres = [10 8];
background = imgaussfilt3(volumeuse, [avgcellradius*5*anisotropyratio]);
dff        = (volumeuse-background)./background;
[dff2, ~]  = spatial_bandpass_3d(volumeuse, avgcellradius, 3, 3, usegpu, anisotropyratio);
background = volumeuse-dff2;
background(background<0) = 0;
%==========================================================================
% detect local maxima and get candidate voxels
smax      = my_max(dff2, 4*sigmause, 1:3);
imgidx    = (dff2==smax) & (dff2 > sdthres(1)); 
imgidx    = gather(imgidx);
imgidx    = imdilate(imgidx, seuse) & gather(dff2 > sdthres(2));

cc    = bwconncomp(imgidx, 18);
cinfo = regionprops3(cc, gather(volumeuse),...
    'PrincipalAxisLength', 'EquivDiameter', 'VoxelList', ...
    'WeightedCentroid','MeanIntensity', 'Solidity');
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
    % %%
    % [~, bb]= pca(zscore([elips,celldiam,cinfo.MeanIntensity]));
    % clf; hold on;
    % [idxcells, C] = kmedoids(bb(:,1:2), 4, 'Replicates',5);
    % for ii =  1:4
    %     plot(bb(idxcells==ii,1),bb(idxcells==ii,2),'o')
    % end
    %%
    ikeep  = ikeepz & ikeepx & ikeepy & ~ilong & ~ismall;
    cinfo  = cinfo(ikeep, :);
    if saveimages
        cim = getCellImages2D(dff, cinfo, sigmause*6);
    end
end
%=========================================================================
% intermediate plotting
if size(cinfo, 1) > 0
    maxc =  sdthres(1);
    
    allvoxels =  cat(1,cinfo.VoxelList{:});
    indtest   = sub2ind(size(dff),allvoxels(:,2), allvoxels(:,1), allvoxels(:,3));
    imgidx = false(size(imgidx));
    imgidx(indtest) = true;
    imshowpair(uint8(max(dff,[],3)*255/maxc),max(imgidx,[],3));
end
%=========================================================================   
end

