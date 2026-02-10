function cell_locations = preprocessColmVolumeBatches(opts)
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

batchsize = 14;
buffsize  = 3;
Nbatches   = ceil(Nslices/batchsize);
snrthres   = 0.8;
normfac = opts.maxdff/(2^16-1);
fid = fopen(opts.fproc,'r');
rangeuse = -15:15;
i0 = 0;
cell_locations = zeros(1e6, 5, 'single');

msg = []; tic;

dx = -10:10;
dy = -10:10;
dz = -2:2;
[xx,yy,zz] = meshgrid(dx, dy, dz);

maxlook = ceil(opts.celldiam./opts.pxsize/2);
cubsize = 2*floor(opts.celldiam./opts.pxsize/2) + 1;
cubsize = [5 5 5];
seuse = strel('cuboid', cubsize);
diamfac = 3*prod(opts.pxsize)/4/pi;

for ibatch = 1:Nbatches
    istart = (ibatch - 1)*batchsize + 1;
    iend   = min(istart + batchsize - 1, Nslices);

    iload  = istart-buffsize:iend+buffsize;
    iload(iload < 1 | iload > Nslices) = [];

    fseek(fid, (iload(1)-1)*Ny*Nx*2, 'bof'); % fseek to batch start in raw file
    dat = fread(fid, [Ny Nx*numel(iload)], '*uint16');
    dat = gpuArray(dat);
    dat = single(dat) * normfac;
    dat = reshape(dat, [Ny, Nx, numel(iload)]);

    % tic;
    % blobs = blobdetect(dat(:,:,1), 7, 'MedianFilter', false, 'GPU', true);
    % toc;


    imgidx = nan(size(dat), 'single');
    for ii = 1:numel(iload)
        smin           = -my_min(-dat(:,:,ii), maxlook([1 2]));
        imgidx(:,:,ii) = (dat(:,:,ii)>smin-1e-3) & (dat(:,:,ii) > snrthres);
    end
    imgidx = imdilate(imgidx, seuse) & gather(dat > (snrthres*0.8));
    cinfo  = bwconncomp(imgidx, 26);

    % smin   = -my_min(-dat, maxlook);
    % imgidx = (dat>smin-1e-3) & (dat > snrthres);
    % imgidx = imdilate(gather(imgidx), seuse) & dat > snrthres/2;
    % cinfo    = bwconncomp(gather(imgidx), 26);

    pxlist   = cinfo.PixelIdxList;
    pxcounts = cellfun(@numel,pxlist);
    ibig     = pxcounts>prod(cubsize)*10;
    ismall   = pxcounts<prod(cubsize)/10;
    pxlist(ibig | ismall) = [];
    ccents = nan(numel(pxlist), 5);
    for ii = 1:numel(pxlist)
        [yidx, xidx, zidx] = ind2sub([Ny Nx numel(iload)], pxlist{ii});
        wts   = dat(pxlist{ii});
        cint  = mean(wts); % get cell intensity
        wts   = wts/sum(wts);
        y0    = yidx' * wts;
        x0    = xidx' * wts;
        z0    = iload(zidx) * wts;
        cdiam = (diamfac*numel(pxlist{ii}))^(1/3); % get cell radius
        ccents(ii, :) = [x0, y0, z0, cdiam, cint];
    end
    % only keep cells within buffer size and far from the edges
    ceval  =  ccents(:, 3) -(iload(1)-1);
    ikeepz = ceval> buffsize & ceval < buffsize+batchsize;
    ikeepx = ccents(:, 1) > Nx * 0.01 & ccents(:, 1) < Nx * 0.99;
    ikeepy = ccents(:, 2) > Ny * 0.01 & ccents(:, 2) < Ny * 0.99;
    ikeep  = ikeepz&ikeepx&ikeepy;
    
    ccents = ccents(ikeep, :);
    % allcells = NaN(size(ccents,1), numel(xx), 'single');
    % 
    % for icell = 1:size(ccents,1)
    %     xcell = round(ccents(icell,1)) + xx;
    %     ycell = round(ccents(icell,2)) + yy;
    %     zcell = round(ccents(icell,3)) + zz -(iload(1)-1);
    %     inds  = sub2ind(size(dat), ycell(:), xcell(:), zcell(:));
    %     volfit = dat(inds);
    %     allcells(icell, :) = volfit;
    %     % tic;
    %     % params = gaussfitn([xcell(:), ycell(:), zcell(:)], double(volfit));
    %     % toc;
    %     % volfit = reshape(volfit, size(xcell));
    %     % imagesc(dx+ccents(icell,1), dy + ccents(icell,2), max(volfit,[],3))
    %     % line(params{3}(1),params{3}(2),'Color','r','Marker','o')
    %     % pause
    % end

    if i0+nnz(ikeep)>size(cell_locations,1)
        cell_locations(1e6 + size(cell_locations,1), 1) = 0;
    end
    cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
    i0 = i0 + nnz(ikeep);

    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Batch %d/%d. Points %d. Time elapsed %2.2f s...\n', ibatch, Nbatches,i0,toc);
    fprintf(msg);
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