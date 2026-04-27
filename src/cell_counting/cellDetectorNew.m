function [cinfo, cim, varargout] = cellDetectorNew(volumeuse, avgcellradius, ...
    sigmause, anisotropyratio, sdthres, bufferzone, saveimages)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
% filterbyprops = getOr(dopts, 'filtercells', true);
filterbyprops = true;
%--------------------------------------------------------------------------
[Ny, Nx, Nz] = size(volumeuse);
usegpu       = isgpuarray(volumeuse);
%--------------------------------------------------------------------------
cubsize = ceil(4*sigmause.*anisotropyratio) + 1;
seuse   = strel('cuboid', cubsize);
%==========================================================================
% perform 3d bandpass filtering
[dff2, ~]  = spatial_bandpass_3d(volumeuse, avgcellradius, 5, 2, usegpu, anisotropyratio);

% get a good estimate of the backgroud
%%
background  = volumeuse-dff2;
isamp       = randperm(numel(background), min(numel(background), 3e4));
currsamp    = background(isamp);
baseline_bg = quantile(currsamp(currsamp > 0), 0.1);
background(background<baseline_bg) = baseline_bg;
backmax = quantile(currsamp(currsamp > 0), 0.9);
background(background>backmax) = backmax;
background  = imgaussfilt3(background, anisotropyratio*2);
% imagesc(max(background,[],3))
%%
%==========================================================================
% the signal is relative to the background
ampsignal = dff2 ./ background;
% imagesc(max(ampsignal,[],3),[0.5 1])
%==========================================================================
% detect local maxima and get candidate voxels
smax      = my_max(ampsignal, 4*sigmause, 1:3);
% candidates should have at least a reasonable intensity
minval    = quantile(dff2(isamp), 0.98, 'all');

imgidx    = (ampsignal==smax) & (ampsignal > sdthres(1)) & (dff2 > minval); 
imgidx    = gather(imgidx);
imgidx    = imdilate(imgidx, seuse) & gather(ampsignal > sdthres(2));

cc    = bwconncomp(imgidx, 18);
cinfo = regionprops3(cc, gather(volumeuse),...
    'PrincipalAxisLength', 'EquivDiameter', 'VoxelList', ...
    'WeightedCentroid','MeanIntensity');
cinfoori = cinfo;
%==========================================================================
cim = [];
if size(cinfo, 1) > 0

    ccents = cinfo.WeightedCentroid;
    %----------------------------------------------------------------------
    % decide which cells to keep
    ikeepz = ccents(:, 3)> bufferzone(3,1) & (ccents(:, 3) < (Nz - bufferzone(3,2)));
    ikeepx = ccents(:, 1)> bufferzone(1,1) & (ccents(:, 1) < (Nx - bufferzone(1,2)));
    ikeepy = ccents(:, 2)> bufferzone(2,1) & (ccents(:, 2) < (Ny - bufferzone(2,2)));
    
    if filterbyprops
        elips    = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
        celldiam = cinfo.EquivDiameter * prod(1./anisotropyratio)^(1/3);
        ilong    = elips>2.5;
        ismall   = celldiam<avgcellradius;
        ibig     = celldiam>avgcellradius*8;
        intthres = 0;
        ihigh    = cinfo.MeanIntensity > intthres;
        ifilter  =  ~ilong & ~ismall & ihigh & ~ibig;
    else
        ifilter = true(size(ccents,1), 1);
    end
    ikeep    = ikeepz & ikeepx & ikeepy & ifilter;
    cinfo    = cinfo(ikeep, :);
    cinfoori = cinfo(ifilter, :);
    %----------------------------------------------------------------------
    if saveimages
        cim = getCellImages2D(dff2, cinfo, sigmause*6);
    end
end

if nargout > 2
    varargout{1} = ampsignal;
end

if nargout > 3
    imgout       = false(Ny, Nx);
    if size(cinfoori, 1) > 0
        allvoxels = cat(1, cinfoori.VoxelList{:});
        indtest   = sub2ind(size(dff2),allvoxels(:,2), allvoxels(:,1), allvoxels(:,3));
        imgidx    = false(size(imgidx));
        imgidx(indtest) = true;
        imgout = max(imgidx,[],3);
    end
    varargout{2} = imgout;
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

