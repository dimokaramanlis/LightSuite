function [backvol, folderproc] = preprocessColmVolumeSimple(dpim, savepath)
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
%%
backvol  = nan(round(scaledown*Ny), round(scaledown*Nx), Nslices, 'single');

maxc = 10; % maximum dff for saving to uint8

Nfilt    = [31 31];
msg = []; tic;
for islice = 1:Nslices

    currim     = imread(fullfile(tim(islice).folder, tim(islice).name));
    backsignal = medfilt2(currim, Nfilt);
    currim     = single(currim);
    backsignal = single(backsignal);
    %------------------------------------------------------------------
    % process and save cell image
    cellimage  = (currim - backsignal)./backsignal;

    % figure;
    % subplot(1,3,1)
    % imagesc(currim,[0 1e3]); axis equal; axis tight;
    % xlim([1500 2200]); ylim([7700 8600])
    % title(islice)
    % colormap(gray)
    % subplot(1,3,2)
    % imagesc(backsignal,[0 1e3]); axis equal; axis tight;
    % xlim([1500 2200]); ylim([7700 8600])
    % title('median filter 21 px')
    % colormap(gray)
    % subplot(1,3,3)
    % imagesc(cellimage,[0 0.5]); axis equal;axis tight;
    % xlim([1500 2200]); ylim([7700 8600])
    % colormap(gray)
    % title('(image - median)/median')


    cellimage(isinf(cellimage) | cellimage<0) = 0;
    cellimage = uint16(((2^16-1)/maxc)*cellimage);
    fname = fullfile(savepath, bpath, sprintf('dff_slice_%04d.tif', indscrr(Nside + 1)));
    imwrite(cellimage, fname,"tif","Compression","lzw")
    %------------------------------------------------------------------
    % for determing edges (mask)
    % 
    %------------------------------------------------------------------
    % fix uniformity of background
    % for background imadjust(imflatfield(backsignal, 1000))
    % findchangepts(std(backvol(:,:,500),[],1),'MaxNumChanges',2,'Statistic','rms')
    % save background
    backvol(:, :, islice) = imresize(backsignal, scaledown);
    
    %------------------------------------------------------------------

    if mod(islice, 2) == 1 || islice == Nslices
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', islice, Nslices,toc);
        fprintf(msg);
    end
end

end