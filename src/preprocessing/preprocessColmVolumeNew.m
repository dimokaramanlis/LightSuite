function [voldown, opts] = preprocessColmVolumeNew(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
tiffpaths   = dir(fullfile(opts.datafolder, '*.tif'));
%--------------------------------------------------------------------------
% other info
scaledownxy  = opts.pxsize(1)/opts.registres;
scaledownz   = opts.pxsize(3)/opts.registres;
[Ny, Nx, Nz] = deal(opts.Ny, opts.Nx,opts.Nz);
medwithfull  = 2*ceil(opts.registres/opts.pxsize(1)/2) + 1;
%--------------------------------------------------------------------------
% initialize collections
backvol  = nan(ceil(scaledownxy*Ny), ceil(scaledownxy*Nx), Nz, 'single');
matuse   = ones(medwithfull);
Nmed     = floor(sum(matuse, 'all')/2);
%--------------------------------------------------------------------------
fid = fopen(opts.fproc, 'W');
msg = []; proctic = tic;


for islice = 1:Nz
    %----------------------------------------------------------------------
    currim = imread(fullfile(tiffpaths(islice).folder, tiffpaths(islice).name));
    % this step is time-limiting 
    filtim = ordfilt2(currim, Nmed, matuse, 'symmetric');
    %----------------------------------------------------------------------
    % write data
    backvol(:, :, islice) = imresize(filtim, scaledownxy, Antialiasing=false);
    % we write a median-filtered version to remove artifacts
    fwrite(fid, medfilt2(currim, [3 3], 'symmetric'), "uint16");
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d. Time per slice %2.2f s. Time elapsed %2.2f s...\n',...
        islice, Nz, toc(proctic)/islice, toc(proctic));
    fprintf(msg);
    %----------------------------------------------------------------------
end
fclose(fid);
%--------------------------------------------------------------------------
fprintf('Saving background volume... '); savetic = tic;

% rescale final volume to match atlas size
voldown             = imresize3(backvol, 'Scale', [1 1 scaledownz]);
voldown(voldown <0) = 0;


% find global threshold
Nxdown    = size(backvol, 2);
T = nan(1,2);
for ii = 1:2
    xvals   = (ii-1)*floor(Nxdown/2) + (1:floor(Nxdown/2));
    currvol = voldown(:, xvals, :);
    maxc    = quantile(voldown(:, xvals, :), 0.999, 'all');
    T(ii)   = maxc*graythresh(currvol/maxc);
end
Tglobal = max(min(T)/2, min(voldown(voldown>0),[],'all'));

% prepare for saving
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
opts.Tglobal    = Tglobal;
% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')


% fpath      = fileparts(opts.fproc);
% regvolpath = fullfile(fpath, sprintf('%s_temporary_reg_volume.dat', opts.mousename));
% 
% fidreg = fopen(regvolpath, 'W');
% fwrite(fidreg, voldown,"uint16");
% fclose(fidreg);

% save(regvolpath, 'voldown')


fprintf('Done! Factor save is %2.1f. Took %2.2f s.\n', regvolfac, toc(savetic));
%--------------------------------------------------------------------------

end