function [finvol, opts] = readSpinalCordSample(dp, sampleres)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

tiffiles = dir(fullfile(dp, '*.tiff'));
tifiles  = dir(fullfile(dp, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
tic;
if numel(tfiles) == 1
    data   = bfopen(fullfile(tfiles.folder, tfiles.name));
    finvol = cat(3, data{1}{:,1});
    Nchans = str2double(data{1}{1,2}(end)); 
    Nz     = size(finvol, 3)/Nchans;
    finvol = reshape(finvol, [size(finvol, [1 2]) Nchans Nz]);
    finvol = permute(finvol, [1 2 4 3]);
else
    Nchans = numel(tfiles);
    finvol = cell(Nchans, 1);
    for ichan = 1 : Nchans
        data          = bfopen(fullfile(tfiles(ichan).folder, tfiles(ichan).name));
        currvol       = cat(3, data{1}{:,1});
        finvol{ichan} = currvol;
    end
    finvol = cat(4, finvol{:});
end
[Nslices, Ny, Nx, Nchan] = size(finvol);
fprintf('Parsed spinal cord sample in %2.1f s. Size %d x %d x %d with %d channels\n', ...
    toc, Nslices, Ny, Nx, Nchan)
%--------------------------------------------------------------------------
opts.datafolder = tfiles(1).folder;
opts.lsfolder   = fullfile(opts.datafolder, 'lightsuite');
makeNewDir(opts.lsfolder)
opts.orisize    = [Nslices, Ny, Nx];
opts.Nchan      = Nchan;
opts.sampleres  = sampleres;
%--------------------------------------------------------------------------



% dataim  = BioformatsImage(dp);
% Nchans = dataim.sizeC;

% finvol  = zeros(dataim.height, dataim.width, dataim.sizeZ, Nchans, 'uint16');
% for ii = 1:Nchans
%     currchan = zeros(dataim.height, dataim.width, dataim.sizeZ, 'uint16');
%     for iz = 1:dataim.sizeZ
%         currchan(:, :, iz) = dataim.getPlane(iz, ii, 1, 1);
%     end
%     currchan(currchan==0) = mode(currchan(currchan>0));
%     finvol(:, :, :, ii) = currchan;
% end

end