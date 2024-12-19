function [volnew, coordsfin] = downsampleLightsheetVolume(dpim, dpsave, pxsize, pxfin)
%DOWNSAMPLELIGHTSHEETVOLUME Summary of this function goes here
%   Detailed explanation goes here

%===========================================================================
makeNewDir(dpsave);
%===========================================================================
bfim = BioformatsImage(dpim);
Nz   = bfim.sizeZ;
Nx   = bfim.width;
Ny   = bfim.height;
Vdata = zeros(Ny, Nx, Nz, 'uint16');
msg = []; tic;
for ii = 1:Nz
    Vdata(:, :, ii) = bfim.getPlane(ii,1,1,1);
    if mod(ii, 20)==1 || ii ==Nz
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', ii, Nz,toc);
        fprintf(msg);
    end
end
%===========================================================================
%%
if size(pxfin) == 1
    pxfin = pxfin * [1 1 1];
end
% downsample volume
voldown      = imresize3(Vdata, 'cubic','Scale', pxsize./pxfin);
[Dy, Dx, Dz] = size(voldown);

finscale = pxfin./pxsize;

xd = 0.5 + finscale(1)*(0:Dx-1)+ finscale(1)/2;
yd = 0.5 + finscale(2)*(0:Dy-1)+ finscale(2)/2;
zd = 0.5 + finscale(3)*(0:Dz-1)+ finscale(3)/2;

coordsfin = {yd(:), xd(:), zd(:)};
save(fullfile(dpsave, 'downcoords.mat'), 'coordsfin')
%===========================================================================
T      = uint16(graythresh(voldown) * single(max(voldown,[],'all')));
maxval = quantile(voldown(voldown>T), 0.99,'all');

volnew = rescale(voldown,"InputMin",T,"InputMax",maxval);
volnew = uint8(rescale(volnew, 0, 255));
%===========================================================================
% save Tiff
[~, fname] = fileparts(dpim);
savepath   = fullfile(dpsave, [fname '_downsampled.tif']);
saveLightsheetVolume(volnew, savepath)
%===========================================================================
end

