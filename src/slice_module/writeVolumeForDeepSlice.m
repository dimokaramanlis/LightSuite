function writeVolumeForDeepSlice(sliceinfo, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

chansmap    = [3 2 1];
namesuse    = sliceinfo.channames;
ncat        = 3-numel(namesuse);
namesuse    = cat(2,  namesuse,repmat(' ',[1 ncat]));

if nargin < 2
    idchan = find(chansmap == 1);
else
    coluse = varargin{1};
    idchan = find(contains(namesuse, coluse));
    idchan = find(idchan == chansmap);
end

errormessage = sprintf('%s ', namesuse{chansmap});
errormessage = sprintf('Select a single channel! Available channels are: %s', errormessage);
assert(~isempty(idchan), errormessage)

fprintf('Channel %s selected to plot for deepslice... ', namesuse{chansmap(idchan)})

dsdir = fullfile(sliceinfo.procpath, 'deepslice');
makeNewDir(dsdir);

volsave = readDownStack(fullfile(sliceinfo.procpath,"volume_for_inspection.tiff"), idchan);

for ii = 1:size(volsave, 3)
    fsavename = sprintf('slice_%s_s%03d.png', sliceinfo.channames{chansmap(idchan)}, ii);
    imwrite(volsave(:,:,ii), fullfile(dsdir, fsavename));
end

fprintf('Done!\n')

end