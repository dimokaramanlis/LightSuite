function [backvol, folderproc] = preprocessColmVolumeTotal(dpim, savepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[~, bpath] = fileparts(dpim);
folderproc = fullfile(savepath, bpath);
makeNewDir(folderproc)
tim     = dir(fullfile(dpim, '*.tif'));
Nslices  = numel(tim);
tfile    = imread(fullfile(tim(1).folder, tim(1).name));
[Ny, Nx] = size(tfile);
Nside    = 2;

pxsize    = [1.44 1.44 5];
scaledown = pxsize(1)/10;
dslice    = -Nside:Nside;
wtsuse    = gausslinefun([0 1 1], dslice);
wtsuse    = wtsuse/sum(wtsuse);
%%
Nbatchx  = 3;
Nbatchy  = 4;
indsx = round(linspace(1, Nx, Nbatchx + 1));
indsy = round(linspace(1, Ny, Nbatchy + 1));

overlap  = 25;
for ix = 1:Nbatchx
    rx = indsx(ix)-overlap:indsx(ix+1)+overlap;
    rx(rx<1 | rx > Nx) = [];
    for iy = 1:Nbatchy
        ry = indsy(iy)-overlap:indsy(iy+1)+overlap;
        ry(ry<1 | ry > Ny) = [];
        currvol = zeros(numel(ry), numel(rx), Nslices, 'uint16');
        for islice = 1:Nslices
            currim = imread(fullfile(tim(islice).folder, tim(islice).name));
            currvol(:,:, islice) = currim(ry, rx);
        end
    end
end




voluse   = nan(Ny, Nx, Nside+1, 'single');
meduse   = nan(Ny, Nx, 2*Nside+1, 'single');
backvol  = nan(round(scaledown*Ny), round(scaledown*Nx), Nslices, 'single');

maxc = 8; % maximum dff for saving to uint8

Nfilt    = [21 21];
msg = []; tic;
for islice = 1:Nslices
    voluse = circshift(voluse, -1, 3);
    meduse = circshift(meduse, -1, 3);
    % meduse(:, :, end) = nan;
    if islice <= Nslices
        currim = imread(fullfile(tim(islice).folder, tim(islice).name));
        filtim = medfilt2(currim, Nfilt);
        voluse(:, :, end) = currim;
        meduse(:, :, end) = filtim;
    end

    if islice >= Nside + 1
        indscrr = (islice - Nside)+ (-Nside:Nside);
        iuse    = indscrr > 0 & indscrr<=Nslices;
        if ~all(iuse)
            wtscurr = wtsuse(iuse);
            wtscurr = wtscurr/sum(wtscurr);
            backsignal = reshape(meduse(:,:,iuse), Ny*Nx, nnz(iuse)) * wtscurr;
        else
            backsignal = reshape(meduse, Ny*Nx, 2*Nside + 1) * wtsuse;
        end
        backsignal = reshape(backsignal, [Ny, Nx]);
        %------------------------------------------------------------------
        % process and save cell image
        cellimage  = (voluse(:, :, 1) - backsignal)./backsignal;

        % cellimage(isinf(cellimage) | cellimage<0) = 0;
        % cellimage = uint8((255/maxc)*cellimage);
        % fname = fullfile(savepath, bpath, sprintf('dff_slice_%04d.tif', indscrr(Nside + 1)));
        % imwrite(cellimage, fname,"tif","Compression","lzw")
        %------------------------------------------------------------------
        % for determing edges (mask)
        % 
        %------------------------------------------------------------------
        % fix uniformity of background
        % for background imadjust(imflatfield(backsignal, 1000))
        % findchangepts(std(backvol(:,:,500),[],1),'MaxNumChanges',2,'Statistic','rms')
        % save background
        backvol(:, :, indscrr(Nside + 1)) = imresize(backsignal, scaledown);
        
        %------------------------------------------------------------------
    end

    
    if mod(islice, 2) == 1 || islice == Nslices
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', islice, Nslices,toc);
        fprintf(msg);
    end
end

end