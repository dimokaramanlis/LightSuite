function cell_locations = extractCellsFromVolume(datavolume, opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%------------------------------------------------------------------------
if isnumeric(datavolume)
    ispath = false;
else
    % it is a path and we have to load it as a path
    assert(isstring(datavolume) | ischar(datavolume))
    ispath = True;
end
%------------------------------------------------------------------------
Ny      = opts.Ny;
Nx      = opts.Nx;
Nslices = opts.Nz;

%------------------------------------------------------------------------
% extract options
normfac    = opts.maxdff/(2^16-1);
sigmause   = [3 3 2];
thresuse   = single([0.5 0.3]);
cellradius = 6;
anisotropy = [1 1 0.6];
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
if ispath
    fid = fopen(datavolume,'r');
end

i0 = 0; itrack = 0; % counters
cell_locations = nan(1e6, 5, 'single');
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
    iloadz     = (istartz-startbuffz):(iendz+endbuffz);


    % we here load data in RAM
    if ispath
        fseek(fid, (iloadz(1)-1)*Ny*Nx*2, 'bof'); % fseek to batch start in raw file
        dat = fread(fid, [Ny Nx*numel(iloadz)], '*uint16');
        dat = reshape(dat, [Ny, Nx, numel(iloadz)]);
    else
        dat = datavolume(:, :, iloadz);
    end
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
            datgpu = single(datgpu) * normfac;
            %----------------------------------------------------------------------
            [ccents, pdiscard] = cellDetector(datgpu, cellradius, sigmause, anisotropy, thresuse);
            if ~isempty(ccents)
                % RECONSIDER THIS, MAYBE WE ARE LOSING SOME CELLS IN EDGES
                ikeepz = ccents(:, 3)> startbuffz & (ccents(:, 3) < (numel(iloadz) - endbuffz));
                ikeepx = ccents(:, 1)> startbuffx & (ccents(:, 1) < (numel(iloadx) - endbuffx));
                ikeepy = ccents(:, 2)> startbuffy & (ccents(:, 2) < (numel(iloady) - endbuffy));
                ikeep  = ikeepz&ikeepx&ikeepy;
                ccents = ccents(ikeep, :);
            
                ccents(:,3) = interp1(1:numel(iloadz), iloadz, ccents(:,3));
                ccents(:,2) = interp1(1:numel(iloady), iloady, ccents(:,2));
                ccents(:,1) = interp1(1:numel(iloadx), iloadx, ccents(:,1));
                
                if i0+nnz(ikeepz)>size(cell_locations,1)
                    cell_locations(1e6 + size(cell_locations,1), 1) = 0;
                end
            
            
                cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
                i0 = i0 + nnz(ikeep);
            end
            %----------------------------------------------------------------------
            % % only keep cells within buffer size and far from the edges
            % cevalz  =  ccents(:, 3) -(iloadz(1)-1);
            % cevalx  =  ccents(:, 1) -(iloadx(1)-1);
            % cevaly  =  ccents(:, 2) -(iloady(1)-1);
            % ikeepz = cevalz> buffsizez & cevalz < buffsizez+batchsizez;
            % ikeepx = cevalx> buffsizexy & cevalx < buffsizexy+batchsizexy;
            % ikeepy = cevaly> buffsizexy & cevaly < buffsizexy+batchsizexy;
            % ikeep  = ikeepz&ikeepx&ikeepy;
            % ccents = ccents(ikeep, :);
            % 
            % if i0+nnz(ikeep)>size(cell_locations,1)
            %     cell_locations(1e6 + size(cell_locations,1), 1) = 0;
            % end
            % cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
            % i0 = i0 + nnz(ikeep);
        
            %----------------------------------------------------------------------
            itrack = itrack + 1;
            fprintf(repmat('\b', 1, numel(msg)));
            msg = sprintf('Batch %d/%d. Points %d. Pdiscard = %2.2f. Time elapsed %2.2f s...\n',...
                itrack, Nbatches,i0,pdiscard*100,toc);
            fprintf(msg);
            %----------------------------------------------------------------------
        end
    end
end
%--------------------------------------------------------------------------
if ispath
    fclose(fid);
end
cell_locations = cell_locations(1:i0, :);
% we finally remove weird entries
irem = any(isnan(cell_locations) | isinf(cell_locations), 2);
cell_locations(irem, :) = [];
%--------------------------------------------------------------------------
% after we are done, save cells
if isfield(opts, 'savepath')
    if ~isempty(opts.savepath)
        makeNewDir(opts.savepath)
        fsavename = fullfile(opts.savepath, 'cell_locations_sample.mat');
        save(fsavename, 'cell_locations')
    end
end
% delete(opts.fproc)

%--------------------------------------------------------------------------
end