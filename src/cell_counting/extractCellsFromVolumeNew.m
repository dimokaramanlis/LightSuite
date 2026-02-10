function cell_locations = extractCellsFromVolumeNew(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%------------------------------------------------------------------------
if opts.debug
    folderdebug = fullfile(opts.savepath, 'cell_detections');
    makeNewDir(folderdebug)
end
%------------------------------------------------------------------------
Ny      = opts.Ny;
Nx      = opts.Nx;
Nslices = opts.Nz;
%------------------------------------------------------------------------
% extract options
sigmause   = [3 3 2];
thresuse   = single([10 8]);
cellradius = 6;
anisotropy = min(opts.pxsize)./opts.pxsize;
voxelvolume = prod(opts.pxsize);
%------------------------------------------------------------------------
% let's figure out batches. TODO: make dependent on cell diameter
batchsizez  = 32;
buffsizez   = ceil(cellradius);
batchsizexy = 1800;
buffsizexy  = ceil(cellradius * 4);

NbatchesZ   = ceil(Nslices/batchsizez);
NbatchesX   = ceil(Nx/batchsizexy);
NbatchesY   = ceil(Ny/batchsizexy);
Nbatches    = NbatchesZ * NbatchesX * NbatchesY;
%------------------------------------------------------------------------
fid = fopen(opts.fproc,'r');

i0 = 0; itrack = 0; % counters
cell_locations = nan(1e6, 6, 'single');
nsigma         = sum(prod(2*6*sigmause(nchoosek(1:3,2))+1,2));

if opts.savecellimages
    cell_images    = nan(1e6, nsigma, 'single');
end
% pback = 0.1;
% if isfield(opts, 'Tglobal')
%     Tglobal = opts.Tglobal;
% else
%     Tglobal = 1;
% end
msg = []; tic;

for ibatchz = 1:NbatchesZ
    %---------------------------------------------------------------------- 
    istartz    = (ibatchz - 1)*batchsizez + 1;
    iendz       = min(istartz + batchsizez - 1, Nslices);
    startbuffz = (ibatchz > 1) * buffsizez;
    endbuffz   = (ibatchz < NbatchesZ) * buffsizez;
    if endbuffz > 0
        endbuffz = min(endbuffz, Nslices - ibatchz*batchsizez);
    end
    iloadz = (istartz-startbuffz):(iendz+endbuffz);

    % dat    = batchLoadSlices(tiffpaths, iloadz, [Ny Nx]);
    fseek(fid, (iloadz(1)-1)*Ny*Nx*2, 'bof'); % fseek to batch start in raw file
    dat = fread(fid, [Ny Nx*numel(iloadz)], '*uint16');
    dat = reshape(dat, [Ny, Nx, numel(iloadz)]);
    %----------------------------------------------------------------------
    % we extract a global threshold
    % ystart    = [1:floor(pback*Ny), Ny-floor(Ny*pback), 1:floor(pback*Ny), Ny-floor(Ny*pback)];
    % xstart    = [1:floor(pback*Nx), Nx-floor(Nx*pback), 1:floor(pback*Nx), Nx-floor(Nx*pback)];
    % valscheck = dat(ystart, xstart, :);
    %----------------------------------------------------------------------
    % we go through X and Y batches if needed
    for ibatchy = 1:NbatchesY
        %----------------------------------------------------------------------
        istarty    = (ibatchy - 1)*batchsizexy + 1;
        iendy      = min(istarty + batchsizexy - 1, Ny);
        startbuffy = (ibatchy > 1) * buffsizexy;
        endbuffy   = (ibatchy < NbatchesY) * buffsizexy;
        if endbuffy > 0
            endbuffy = min(endbuffy, Ny - ibatchy*batchsizexy);
        end
        iloady     = (istarty-startbuffy):(iendy+endbuffy);
        %----------------------------------------------------------------------
        for ibatchx = 1:NbatchesX
            %----------------------------------------------------------------------
            istartx    = (ibatchx - 1)*batchsizexy + 1;
            iendx      = min(istartx + batchsizexy - 1, Nx);
            startbuffx = (ibatchx > 1) * buffsizexy;
            endbuffx   = (ibatchx < NbatchesX) * buffsizexy;
            if endbuffx > 0
                endbuffx = min(endbuffx, Nx - ibatchx*batchsizexy);
            end
            iloadx     = (istartx-startbuffx):(iendx+endbuffx);
            %----------------------------------------------------------------------
            % data goes into the gpu
            datgpu = gpuArray(dat(iloady, iloadx, :));
            datgpu = single(datgpu);
            %----------------------------------------------------------------------
            % extract candidate cells with info
            bufferzone   = [startbuffx endbuffx; startbuffy endbuffy; startbuffz endbuffz];
            [cinfo, cim] = cellDetectorNew(datgpu, cellradius, ...
                sigmause, anisotropy, thresuse, bufferzone, opts.savecellimages);
            if ~isempty(cinfo)
                ccents = cinfo.WeightedCentroid;
                %----------------------------------------------------------
                ccents(:,3) = interp1(1:numel(iloadz), iloadz, ccents(:,3));
                ccents(:,2) = interp1(1:numel(iloady), iloady, ccents(:,2));
                ccents(:,1) = interp1(1:numel(iloadx), iloadx, ccents(:,1));
                %----------------------------------------------------------
                eqdiam    = cinfo.EquivDiameter * voxelvolume^(1/3);
                intmean   = cinfo.MeanIntensity;

                elratio   = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
                cellfeats = [intmean eqdiam elratio];
                % elsort  = sort(elaxes, 2, 'descend');
                % elinds  = [elsort(:, 2)./elsort(:,1), elsort(:,3)./elsort(:,2)];
                % 
                % 
                % fanum   = sqrt(sum((elaxes - mean(elaxes, 2)).^2, 2));
                % faden   = sqrt(sum(elaxes.^2, 2));
                % elfa    = fanum./faden;
                %----------------------------------------------------------
                if i0+size(ccents, 1)>size(cell_locations,1)
                    cell_locations(1e6 + size(cell_locations,1), 1) = 0;
                    if opts.savecellimages
                        cell_images(1e6 + size(cell_locations,1), 1) = 0;
                    end
                end

                cell_locations(i0 + (1:size(ccents, 1)), :) = [ccents cellfeats];
                
                if opts.savecellimages
                    cell_images(i0 + (1:size(ccents, 1)), :) = cim;
                end
                i0 = i0 + size(ccents, 1);
            end
            %----------------------------------------------------------------------
            if opts.debug & size(ccents,1) > 200
                pathslice = fullfile(folderdebug, ...
                    sprintf('%03d_batch_%d_detections.png', itrack, size(ccents,1)));
                imwrite(imgout, pathslice,"png","BitDepth",8)
            end           
            %----------------------------------------------------------------------
            itrack = itrack + 1;
            fprintf(repmat('\b', 1, numel(msg)));
            msg = sprintf('Batch %d/%d. Points %d. Time elapsed %2.2f s...\n',...
                itrack, Nbatches,i0,toc);
            fprintf(msg);
            %----------------------------------------------------------------------
        end
    end
end
%--------------------------------------------------------------------------
fclose(fid);

cell_locations = cell_locations(1:i0, :);

% we finally remove weird entries
irem = any(isnan(cell_locations) | isinf(cell_locations), 2);
cell_locations(irem, :) = [];
if opts.savecellimages
    cell_images    = cell_images(1:i0, :);
    cell_images(irem, :)    = [];
end

%--------------------------------------------------------------------------
% after we are done, save cells
if isfield(opts, 'savepath')
    if ~isempty(opts.savepath)
        makeNewDir(opts.savepath)
        fsavename = fullfile(opts.savepath, 'cell_locations_sample.mat');
        save(fsavename, 'cell_locations')
        if opts.savecellimages
            imwindow = sigmause*6;
            save(fsavename, 'cell_locations', 'cell_images', 'imwindow')
        end
    end
end
% delete(opts.fproc)

%--------------------------------------------------------------------------
end


% 
% 
% 
% 
%                 X = [elips,cinfo.MeanIntensity,cinfo.Solidity, cinfo.EquivDiameter];
% 
% 
% prem   = 0;
% cellimages = zeros(0, prod(2*4*sigmause(1:2)+1));
% if ~isempty(cinfo)
% 
%     elips   = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
%     ilong   = elips>2.5;
%     ismall  = cinfo.EquivDiameter<cellradius/2;
%     ilow    = cinfo.MeanIntensity < quantile(cinfo.MeanIntensity,0.99)/2;
% 
%     X = [elips,cinfo.MeanIntensity,cinfo.Solidity, cinfo.EquivDiameter];
%     [aa,bb] = pca(zscore(X));
%     Nclust = 6;
%     [idx, C] = kmeans(bb(:,1:3), Nclust, 'Replicates', 10,'MaxIter',200);
%     clf;
%     for ii = 1:size(C,1)
%         hold on;
%         scatter3(bb(idx==ii,1),bb(idx==ii,2),bb(idx==ii,3))
%     end
%     diamsall = accumarray(idx, X(:,4),[],@median);
%     elipsall = accumarray(idx, X(:,1),[],@median);
%     solall   = accumarray(idx,  X(:,3),[],@median);
%     intall   = accumarray(idx,  X(:,2),[],@median);
%     ikeep   = diamsall > max(diamsall)*0.5 & elipsall < 2.5 & solall<0.99...
%         & intall > max(intall)*0.2;
%     iweird  = ~ismember(idx, find(ikeep));
% 
% 
% 
%     % iweird = ilong | ismall | ilow;
%     % prem   = nnz(iweird)/size(cinfo, 1);
%     % 
%     % if nnz(~iweird)>0
%     %     allvoxels =  cat(1,cinfo(~iweird,:).VoxelList{:});
%     %     indtest = sub2ind([Ny, Nx],allvoxels(:,2), allvoxels(:,1));
%     %     imgout(indtest) = true;
%     % 
%     %     imtosave = max(dff, [], 3);
%     %     imtosave = gather(uint8(255 * imtosave/thresSNR(1)));
%     %     imtosave(imgout) = 255;
%     %     imtosave = cat(3, uint8(imgout*255), imtosave, uint8(imgout*255));
%     % end
% 
% 
%     allvoxels =  cat(1,cinfo(~iweird,:).VoxelList{:});
%     indtest = sub2ind(size(volumeuse),allvoxels(:,2), allvoxels(:,1), allvoxels(:,3));
% 
%     imgidx = false(size(imgidx));
%     imgidx(indtest) = true;
%     maxc = quantile(dff2, 0.999, 'all');
%     imshowpair(uint8(max(dff2,[],3)*255/maxc),max(imgidx,[],3));
%     idplot = 15:25; imshowpair(uint8(max(dff2(:,:,idplot),[],3)*255/maxc),max(imgidx(:,:,idplot),[],3));
% 
%     if saveim
%         cellimages = getCellImages2D(volumeuse, cinfo, sigmause*6);
%         cellimages(iweird, :) = [];
%     end
%     cinfo(iweird, :) = [];    
% 
% end
% %--------------------------------------------------------------------------
% % package cell properties
% cpoints    = [cinfo.WeightedCentroid cinfo.EquivDiameter cinfo.MeanIntensity];
% 
% 
% 
% 
% 
% 
