function [finvol, opts] = readSpinalCordSample(dp, sampleres)
%READSPINALCORDSAMPLE Load a multi-channel spinal cord image volume from disk.
%
%   NOTE: the registration pipeline no longer goes through this function. Cord
%   samples are read by READLIGHTSHEETOPTS and downsampled by
%   PREPROCESSLIGHTSHEETVOLUME, exactly like brains, so that both share one
%   loader and one set of options. This is kept for pulling a whole cord volume
%   into memory outside the pipeline, e.g. for inspection.
%
%   [FINVOL, OPTS] = READSPINALCORDSAMPLE(DP, SAMPLERES) reads the TIFF
%   image data stored in folder DP and returns the 4-D image volume FINVOL
%   together with a struct OPTS describing the sample and the intended
%   registration.
%
%   Two on-disk layouts are supported automatically:
%     * A single .tif/.tiff file containing a multi-channel z-stack. The
%       channels are read one at a time with BioformatsImage.
%     * One .tif/.tiff file per channel. Each file is read with BFOPEN and
%       the channels are concatenated along the 4th dimension.
%
%   In both cases, zero-valued voxels (typically background/padding produced
%   during acquisition or stitching) are replaced by the mode of the
%   non-zero voxels in that channel, so empty regions take a representative
%   background intensity instead of a hard zero.
%
%   INPUTS
%     dp        - Path to the folder containing the .tif/.tiff data.
%     sampleres - 1x3 voxel resolution of the acquired sample. Stored
%                 verbatim in OPTS.sampleres; not otherwise interpreted here.
%
%   OUTPUTS
%     finvol - uint16 array of size [Ny, Nx, Nz, Nchan] holding the image
%              volume for every channel.
%     opts   - Struct with fields:
%                datafolder      - Folder the data was read from.
%                lsfolder        - 'lightsuite' subfolder (created if absent).
%                orisize         - [Ny, Nx, Nz] original volume size.
%                Nchan           - Number of channels.
%                sampleres       - Voxel resolution of the input (= SAMPLERES).
%                registrationres - Target voxel resolution for registration.
%
%   See also BIOFORMATSIMAGE, BFOPEN, MAKENEWDIR.

% Collect .tif and .tiff files (tif first, then tiff, as before).
tfiles = [dir(fullfile(dp, '*.tif')); dir(fullfile(dp, '*.tiff'))];

tic;
if isscalar(tfiles)
    % --- Single multi-channel stack -------------------------------------
    dataim = BioformatsImage(fullfile(tfiles.folder, tfiles.name));
    Nchan  = dataim.sizeC;

    finvol = zeros(dataim.height, dataim.width, dataim.sizeZ, Nchan, 'uint16');
    for ic = 1:Nchan
        vol = zeros(dataim.height, dataim.width, dataim.sizeZ, 'uint16');
        for iz = 1:dataim.sizeZ
            vol(:, :, iz) = dataim.getPlane(iz, ic, 1, 1);
        end
        finvol(:, :, :, ic) = fillZerosWithMode(vol);
    end
else
    % --- One file per channel -------------------------------------------
    Nchan = numel(tfiles);
    chans = cell(Nchan, 1);
    for ic = 1:Nchan
        data      = bfopen(fullfile(tfiles(ic).folder, tfiles(ic).name));
        chans{ic} = fillZerosWithMode(cat(3, data{1}{:, 1}));
    end
    finvol = cat(4, chans{:});
end

% Note: dim 1 = rows (Ny), dim 2 = cols (Nx), dim 3 = slices (Nz).
[Ny, Nx, Nz, Nchan] = size(finvol);
fprintf('Parsed spinal cord sample in %2.1f s. Size %d x %d x %d with %d channels\n', ...
    toc, Ny, Nx, Nz, Nchan)

% --- Assemble output options --------------------------------------------
opts.datafolder      = tfiles(1).folder;
opts.lsfolder        = fullfile(opts.datafolder, 'lightsuite');
makeNewDir(opts.lsfolder);
opts.orisize         = [Ny, Nx, Nz];
opts.Nchan           = Nchan;
opts.sampleres       = sampleres;
opts.registrationres = [20 20 20];   % target resolution for registration

end

% =========================================================================
function vol = fillZerosWithMode(vol)
%FILLZEROSWITHMODE Replace zero voxels with the mode of the non-zero voxels.
vol(vol == 0) = mode(vol(vol > 0));
end