function cell_locations = extractCellsFromVolume(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 
% [~, bpath] = fileparts(dpim);
% % folderproc = fullfile(savepath, bpath);
% % makeNewDir(folderproc)
% tim     = dir(fullfile(dpim, '*.tif'));
% Nslices  = numel(tim);
% tfile    = imread(fullfile(tim(1).folder, tim(1).name));
% [Ny, Nx] = size(tfile);
Ny      = opts.Ny;
Nx      = opts.Nx;
Nslices = opts.Nz;

%%

batchsizez  = 32;
buffsizez   = 6;
NbatchesZ   = ceil(Nslices/batchsizez);

batchsizexy = 2000;
buffsizexy  = 24;
NbatchesX   = ceil(Nx/batchsizexy);
NbatchesY   = ceil(Ny/batchsizexy);

snrthres   = 0.8;
normfac = opts.maxdff/(2^16-1);


fid = fopen(opts.fproc,'r');
i0 = 0;
cell_locations = zeros(1e6, 5, 'single');

msg = []; tic;

dx = -10:10;
dy = -10:10;
dz = -2:2;
[xx,yy,zz] = meshgrid(dx, dy, dz);

maxlook = ceil(opts.celldiam./opts.pxsize/2);
cubsize = 2*floor(opts.celldiam./opts.pxsize/2) + 1;

diamfac = 3*prod(opts.pxsize)/4/pi;

sigmause    = ceil(opts.celldiam./opts.pxsize/4);
sigmause(3) = sigmause(3) * 2;
fsuse = fspecial3('log',2*ceil(sigmause*2)+1, sigmause);
fsuse = fsuse./norm(fsuse(:));

cubsize = 2*sigmause ;
seuse = strel('cuboid', cubsize);

sigmasall = [sigmause/2;sigmause;sigmause*2];

Thres    = [2 1];
ThresSNR = [0.6 0.4];
itrack = 0;
Nbatches = NbatchesZ * NbatchesX * NbatchesY;

dty = -10:10;
dtx = -10:10;
dtz = -5:5;

for ibatchz = 1:NbatchesZ
    ibatchz = 25;
    %----------------------------------------------------------------------    
    istart = (ibatchz - 1)*batchsizez + 1;
    iend   = min(istart + batchsizez - 1, Nslices);

    iload  = istart-buffsizez:iend+buffsizez;
    iload(iload < 1 | iload > Nslices) = [];
    % we here load data in RAM
    fseek(fid, (iload(1)-1)*Ny*Nx*2, 'bof'); % fseek to batch start in raw file
    dat = fread(fid, [Ny Nx*numel(iload)], '*uint16');
    dat = reshape(dat, [Ny, Nx, numel(iload)]);
    %----------------------------------------------------------------------
    for ibatchy = 1:NbatchesY
        ibatchy = 4;
        istarty = (ibatchy - 1)*batchsizexy + 1;
        iendy   = min(istarty + batchsizexy - 1, Ny);

        iloady  = istarty-buffsizexy:iendy+buffsizexy;
        iloady(iloady < 1 | iloady > Ny) = [];

        for ibatchx = 1:NbatchesX
            ibatchx = 3;
            istartx = (ibatchx - 1)*batchsizexy + 1;
            iendx   = min(istartx + batchsizexy - 1, Nx);

            iloadx  = istartx-buffsizexy:iendx+buffsizexy;
            iloadx(iloadx < 1 | iloadx > Nx) = [];
            %----------------------------------------------------------------------
            datgpu = gpuArray(dat(iloady, iloadx, :));
            datgpu = single(datgpu) * normfac;
            % first filter for cleaning 
            % THIS FILTERING SHOULD CHANGE

            % 
            % 
            % % let's learn cells from data
            [fout, power_retained_spectrum] = spatial_bandpass_3d(datgpu, 6, ...    
                4, 3, true, [1 1 0.6]);
            % find maxima in filtered space
            smin           = -my_min(-fout, sigmause, 3);
            % we keep maxima with a minimum snr
            % this snr should be probably in filtered space
            imgidx = (fout>smin-1e-3) & (datgpu > ThresSNR(1)) & fout>0;
            imgidx = gather(imgidx);
            imgidx = imdilate(imgidx, seuse) & gather(datgpu > ThresSNR(2));
            cinfo  = regionprops3(imgidx, gather(datgpu),'all');

            alleigs = cat(2,cinfo.EigenValues{:});
            el1     = sqrt(1 - alleigs(3,:)./alleigs(1,:));
            el2     = sqrt(1 - alleigs(3,:)./alleigs(2,:));

            % 
            X = nan(size(cinfo,1), numel(dtx)*numel(dty)*numel(dtz),'single');
            for icell = 1:size(cinfo,1)
                ccent  =  round(cinfo.Centroid(icell,:));
                xind = ccent(1) + dtx;
                yind = ccent(2) + dty;
                zind = ccent(3) + dtz;
                if any(zind<1)|any(xind<1)|any(yind<1)|any(xind>batchsizexy)|any(yind>batchsizexy)|...
                        any(zind>numel(iload))
                    continue
                end
                currvol = fout(yind, xind, zind);
                X(icell,:) = reshape(currvol, 1, []);
            end
            irem = any(isnan(X),2);
            XX = X(~irem, :);
            XX = XX./sqrt(sum(XX.^2,2));
            [aa,bb,cc] = svd(XX);
            aplot = reshape(cc(:,1),numel(dty),numel(dtx),numel(dtz));
            % 
            % % newdata = imgidx.*datgpu;
            % cinfo1  = bwconncomp(imgidx, 26);
            % 
            % aplot = aplot./norm(aplot(:));
            % foutnew = imfilter(datgpu, gpuArray(double(aplot)), 'replicate', 'same');
            % 
            % for ii = 1:3
            %     fscurr = fspecial3('log',3*sigmause, sigmasall(ii,:));
            %     norm(fscurr(:))
            %     fscurr = fscurr./norm(fscurr(:));
            %     foutcurr = imfilter(datgpu, -gpuArray(fscurr), 'replicate', 'same');
            %     imagesc(foutcurr(:,:,11),[0 2])
            %     pause;
            % end
            tic;
            fout = imfilter(datgpu, -gpuArray(fsuse), 'replicate', 'same');
            toc;

            % then we find 3d maxima in filtered space
            smin           = -my_min(-fout, sigmause, 3);
            % we keep maxima with a minimum snr
            % this snr should be probably in filtered space
            imgidx = (fout>smin-1e-3) & (fout > 1) & (datgpu > ThresSNR(1));
            imgidx = gather(imgidx);
            imgidx = imdilate(imgidx, seuse) & gather((datgpu > ThresSNR(2)));
            % newdata = imgidx.*datgpu;
            cinfo  = bwconncomp(imgidx, 26);
            cinfo  = regionprops3(imgidx, gather(datgpu),'all');
            % infouse = regionprops(cinfo);
            
            pxlist   = cinfo.PixelIdxList;
            pxcounts = cellfun(@numel,pxlist);
            ibig     = pxcounts>prod(cubsize)*10;
            ismall   = pxcounts<prod(cubsize)/10;
            pxlist(ibig | ismall) = [];
            ccents = nan(numel(pxlist), 5);
            for ii = 1:numel(pxlist)
                [yidx, xidx, zidx] = ind2sub([numel(iloady) numel(iloadx) numel(iload)],...
                    pxlist{ii});
                wts   = datgpu(pxlist{ii});
                cint  = mean(wts); % get cell intensity
                wts   = wts/sum(wts);
                y0    = iloady(yidx) * wts;
                x0    = iloadx(xidx) * wts;
                z0    = iload(zidx) * wts;
                cdiam = (diamfac*numel(pxlist{ii}))^(1/3); % get cell radius
                ccents(ii, :) = [x0, y0, z0, cdiam, cint];
            end
            % only keep cells within buffer size and far from the edges
            cevalz  =  ccents(:, 3) -(iload(1)-1);
            cevalx  =  ccents(:, 1) -(iloadx(1)-1);
            cevaly  =  ccents(:, 2) -(iloady(1)-1);
            ikeepz = cevalz> buffsizez & cevalz < buffsizez+batchsizez;
            ikeepx = cevalx> buffsizexy & cevalx < buffsizexy+batchsizexy;
            ikeepy = cevaly> buffsizexy & cevaly < buffsizexy+batchsizexy;
            ikeep  = ikeepz&ikeepx&ikeepy;
            ccents = ccents(ikeep, :);

            if i0+nnz(ikeep)>size(cell_locations,1)
                cell_locations(1e6 + size(cell_locations,1), 1) = 0;
            end
            cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
            i0 = i0 + nnz(ikeep);
        
            itrack = itrack + 1;
            fprintf(repmat('\b', 1, numel(msg)));
            msg = sprintf('Batch %d/%d. Points %d. Time elapsed %2.2f s...\n', itrack, Nbatches,i0,toc);
            fprintf(msg);
            % for every maximum, we ask the classifier
            %----------------------------------------------------------------------
        end
    end
end
fclose(fid);
cell_locations = cell_locations(1:i0, :);
%--------------------------------------------------------------------------
% after we are done, save cells and delete binary file
makeNewDir(opts.savepath)
fsavename = fullfile(opts.savepath, 'cell_locations_sample.mat');
save(fsavename, 'cell_locations')
delete(opts.fproc)
%--------------------------------------------------------------------------
end