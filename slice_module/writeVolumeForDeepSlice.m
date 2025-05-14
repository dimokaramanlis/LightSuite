function writeVolumeForDeepSlice(sliceinfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dsdir = fullfile(sliceinfo.procpath, 'deepslice');
makeNewDir(dsdir);

volsave = readDownStack(fullfile(sliceinfo.procpath,"volume_for_inspection.tiff"), 1);

for ii = 1:size(volsave, 3)
    fsavename = sprintf('slice_s%03d.png', ii);
    imwrite(volsave(:,:,ii), fullfile(dsdir, fsavename));
end


end