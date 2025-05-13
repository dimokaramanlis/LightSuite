function saveLargeSliceVolume(volsave, channames, folderpath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nchannels = numel(channames);
makeNewDir(folderpath)

options.compress = 'lzw';
options.message  = false;
options.color    = false;
options.big      = true;

for ichan = 1:Nchannels
    currname = sprintf('chan%02d_%s.tiff',ichan, channames{ichan});
    volpath  = fullfile(folderpath, currname);
    if exist(volpath, 'file')
        delete(volpath);
    end
    saveastiff(squeeze(volsave(:,:,ichan,:)), volpath, options);
end



end