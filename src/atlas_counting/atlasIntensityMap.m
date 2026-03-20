function mousemap = atlasIntensityMap(mousecells, avsize, sigma, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if iscell(mousecells)
    Nmice     = numel(mousecells);
else
    Nmice = 1;
    mousecells = {mousecells};
end

if nargin < 4
    wts = ones(Nmice, 1)/Nmice;
else
    wts = varargin{1}(:);
    wts = wts/sum(wts);
end

Nmid         = avsize(3)/2;
mousemap     = zeros([avsize(1:2) Nmid], 'single');
[Ny,Nx,Nz]   = size(mousemap);
[xx,yy,zz]   = meshgrid([1:100 Nx-100:Nx], 1:Ny, 1:100);
indsback     = sub2ind(size(mousemap), yy(:), xx(:), zz(:));

for ii = 1:Nmice 
    registeredvolume = getIntensityVolume(mousecells{ii});
    valsback         = getRandomValsFromArray(registeredvolume(indsback), 2e4);
    backval          = single(median(valsback(valsback>0)));
    mousemap         = mousemap + ((single(registeredvolume)-backval)/backval) * wts(ii);
end

% avplot       = (avplot(:,:, 1:Nmid) + flip(avplot(:,:, (Nmid + 1):end), 3))/2;
if sigma > 0
    mousemap = imgaussfilt3(mousemap, sigma);
end

end

function regvol = getIntensityVolume(currpath)

savepath         = fullfile(currpath, 'volume_atlas_space.tiff');
if exist(savepath, "file")
    regvol = readDownStack(savepath);
else
    transform_params = load(fullfile(currpath, 'transform_params.mat'));
    backvolume       = readDownStack(fullfile(currpath, 'sample_register_20um.tif'));
    regvol           = atlasSpaceFromVolumeParams(backvolume, transform_params);
    Nmid             = size(regvol, 3)/2;
    regvol           = single(regvol);
    regvol           = (regvol(:,:, 1:Nmid) + flip(regvol(:,:, (Nmid + 1):end), 3))/2;
    regvol           = uint16(regvol);
    saveLightsheetVolume(regvol, savepath, 16);
end




end