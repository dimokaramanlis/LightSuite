function cell_locations = extractCellsFromVolumeNew(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%------------------------------------------------------------------------
opts.prefix = getOr(opts, 'prefix', '');
writetocsv  = getOr(opts, 'writetocsv', false);
usegpu      = getOr(opts, 'usegpu', true);
%------------------------------------------------------------------------
if opts.debug
    fprintf('Using debug mode, you will get pictures with cell detections\n');
    folderdebug = fullfile(opts.savepath, sprintf('%scell_detections', opts.prefix));
    if exist(folderdebug,'dir')
        rmdir(folderdebug, 's')
    end
    makeNewDir(folderdebug);
end
%------------------------------------------------------------------------
Ny      = opts.Ny;
Nx      = opts.Nx;
Nslices = opts.Nz;
%------------------------------------------------------------------------
% extract options
thresuse    = single(getOr(opts, 'thres_cell_detect', [0.5 0.4]));
if thresuse(2) > thresuse(1)
    warning(['The second detection threshold should be smaller or equal than the first.'...
        ' Setting both to %2.2f.'], thresuse(1))
    thresuse(2) = thresuse(1);
end

voxelvolume = prod(opts.pxsize);
cellradius  = ceil(opts.celldiam/(2*min(opts.pxsize)) + 1);
anisotropy  = min(opts.pxsize)./opts.pxsize;
sigmause    = max(ceil(anisotropy.*cellradius/2), 2);
%------------------------------------------------------------------------
% let's figure out batches.
batchsizez  = getOr(opts, 'batchsizez', 32);
batchsizexy = getOr(opts, 'batchsizexy', 1800);
buffsizez   = floor(cellradius*anisotropy(3)*4);
buffsizexy  = floor(cellradius*anisotropy(1)*4);

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

% we also triage cells based on their proximity to edges to get rid of
% cells outside the tissue
mind    = getOr(opts, 'minedgedist', 50);
sizevol = [opts.Nx, opts.Ny, opts.Nz].*opts.pxsize;

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
            if usegpu
                datgpu = gpuArray(dat(iloady, iloadx, :));
            else
                 datgpu = dat(iloady, iloadx, :);
            end
            datgpu = single(datgpu);
            %----------------------------------------------------------------------
            % extract candidate cells with info
            bufferzone   = [startbuffx endbuffx; startbuffy endbuffy; startbuffz endbuffz];
            [cinfo, cim, ampmax, imgout] = cellDetectorNew(datgpu, cellradius, ...
                sigmause, anisotropy, thresuse, bufferzone, opts.savecellimages);
            if ~isempty(cinfo) % we first check locations
                ccents = cinfo.WeightedCentroid;
                %----------------------------------------------------------
                ccents(:,3) = interp1(1:numel(iloadz), iloadz, ccents(:,3));
                ccents(:,2) = interp1(1:numel(iloady), iloady, ccents(:,2));
                ccents(:,1) = interp1(1:numel(iloadx), iloadx, ccents(:,1));
                %----------------------------------------------------------
                % adding triage for volume edges
                locsum  = ccents.*opts.pxsize;
                iedge   = any(locsum< mind,2) | any(locsum > (sizevol-mind), 2);
                ccents  = ccents(~iedge, :);
                cinfo   = cinfo(~iedge,  :);
                if opts.savecellimages
                    cim = cim(~iedge, :);
                end
                %----------------------------------------------------------
            end
            if ~isempty(cinfo) % we then package everything
                %----------------------------------------------------------
                eqdiam    = cinfo.EquivDiameter * voxelvolume^(1/3);
                intmean   = cinfo.MeanIntensity;

                elratio   = cinfo.PrincipalAxisLength(:,1)./cinfo.PrincipalAxisLength(:,2);
                cellfeats = [intmean eqdiam elratio];
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
                %----------------------------------------------------------
                if opts.debug & size(ccents,1) > 1
                    pathslice = fullfile(folderdebug, ...
                        sprintf('%03d_batch_x%03d_y%03d_z%03d_%d_detections.png', ...
                        itrack, ibatchx,ibatchy,ibatchz, size(ccents,1)));
                    ampmax   = max(ampmax,[], 3);
                    imtosave = gather(uint8(255 * ampmax/thresuse(1)));
                    imtosave(imgout) = 255;
                    imtosave = cat(3, uint8(imgout*255), imtosave, uint8(imgout*255));
                    imtosave = imresize(imtosave, 0.5);
                    imwrite(imtosave, pathslice,"png","BitDepth",8)
                end
                %----------------------------------------------------------
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
%%
%--------------------------------------------------------------------------
fclose(fid);
%--------------------------------------------------------------------------
cell_locations = cell_locations(1:i0, :);

% we finally remove weird entries
irem                    = any(isnan(cell_locations) | isinf(cell_locations), 2);
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
        fsavename = fullfile(opts.savepath, sprintf('%scell_locations_sample.mat', opts.prefix));
        save(fsavename, 'cell_locations')
        if opts.savecellimages
            imwindow = sigmause*6;
            save(fsavename, 'cell_locations', 'cell_images', 'imwindow', '-v7.3')
        end
        if writetocsv
            csvfilename = fullfile(opts.savepath, sprintf('%scell_locations_sample.csv', opts.prefix));
            writematrix(cell_locations, csvfilename);
        end
    end
end
%--------------------------------------------------------------------------
end
