function [voldown, opts] = preprocessColmVolume(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
tiffpaths   = dir(fullfile(opts.datafolder, '*.tif'));
%--------------------------------------------------------------------------
% how to combine slices over z
Nsidez       = ceil(opts.celldiam/opts.pxsize(3));
dslice       = -Nsidez:Nsidez;
wtsuse       = gausslinefun([0 1 1], dslice);
wtsuse       = wtsuse/sum(wtsuse);
%--------------------------------------------------------------------------
% other info
scaledownxy  = opts.pxsize(1)/opts.registres;
scaledownz   = opts.pxsize(3)/opts.registres;
[Ny, Nx, Nz] = deal(opts.Ny, opts.Nx,opts.Nz);
diaminpx     = opts.celldiam./opts.pxsize(1);
medwithfull  = 2*ceil((2.5*diaminpx)/2) + 1;
%--------------------------------------------------------------------------
% initialize collections
voluse   = nan(Ny, Nx, Nsidez+1,   'single');
meduse   = nan(Ny, Nx, 2*Nsidez+1, 'single');
backvol  = nan(ceil(scaledownxy*Ny), ceil(scaledownxy*Nx), Nz, 'single');
maxc     = opts.maxdff; % maximum dff for saving to uint8
normfac  = single(intmax("uint16"))/maxc;

matuse   = ones(medwithfull);
Nmed     = floor(sum(matuse, 'all')/2);
%--------------------------------------------------------------------------
fid = fopen(opts.fproc, 'W');
msg = []; tic;


for islice = 1:Nz+Nsidez
    voluse = circshift(voluse, -1, 3);
    meduse = circshift(meduse, -1, 3);
    % meduse(:, :, end) = nan;
    if islice <= Nz
        currim = imread(fullfile(tiffpaths(islice).folder, tiffpaths(islice).name));
        
        % this step is time-limiting 

        filtim            = ordfilt2(currim, Nmed, matuse);
        voluse(:, :, end) = medfilt2(currim, [3 3]); % to remove speckles and line artifacts; 
        meduse(:, :, end) = filtim;
    end

    if islice >= Nsidez + 1
        indscrr = (islice - Nsidez)+ (-Nsidez:Nsidez);
        iuse    = indscrr > 0 & indscrr<=Nz;
        if ~all(iuse)
            wtscurr    = wtsuse(iuse);
            wtscurr    = wtscurr/sum(wtscurr);
            backsignal = reshape(meduse(:,:,iuse), Ny*Nx, nnz(iuse)) * wtscurr;
        else
            backsignal = reshape(meduse, Ny*Nx, 2*Nsidez + 1) * wtsuse;
        end
        backsignal = reshape(backsignal, [Ny, Nx]);
        %------------------------------------------------------------------
        % process and save cell image
        cellimage  = (voluse(:, :, 1) - backsignal)./backsignal;
        cellimage(isinf(cellimage) | cellimage<0) = 0;
        cellimage = uint16(normfac*cellimage);

        fwrite(fid, cellimage, "uint16");

        % fname = fullfile(savepath, bpath, sprintf('dff_slice_%04d.tif', indscrr(Nside + 1)));
        % imwrite(cellimage, fname,"tif","Compression","packbits")
        %------------------------------------------------------------------
        % save background
        backvol(:, :, indscrr(Nsidez + 1)) = imresize(backsignal, scaledownxy);
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d. Time per slice %2.2f s. Time elapsed %2.2f s...\n',...
        islice, Nz, toc/islice, toc);
    fprintf(msg);
    %----------------------------------------------------------------------
end
fclose(fid);
%--------------------------------------------------------------------------
fprintf('Saving background volume... ');tic;

% rescale final volume to match atlas size
voldown   = imresize3(backvol, 'Scale', [1 1 scaledownz]);
regvolfac = (2^16-1)/max(voldown, [],"all");
voldown   = uint16(regvolfac*voldown);


% save volume for control point and registration
samplepath = fullfile(opts.savepath, sprintf('sample_register_%dum.tif', opts.registres));
options.compress = 'lzw';
options.message  = false;
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(voldown, samplepath, options);



opts.regvolpath = samplepath;
opts.regvolfac  = regvolfac;
opts.regvolsize = size(voldown);

% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')


% fpath      = fileparts(opts.fproc);
% regvolpath = fullfile(fpath, sprintf('%s_temporary_reg_volume.dat', opts.mousename));
% 
% fidreg = fopen(regvolpath, 'W');
% fwrite(fidreg, voldown,"uint16");
% fclose(fidreg);

% save(regvolpath, 'voldown')


fprintf('Done! Factor save is %2.1f. Took %2.2f s.\n', regvolfac, toc);
%--------------------------------------------------------------------------
end