    % dp            = 'D:\example_charlie\CGF028'; % OG path, good
% dp             = 'D:\example_charlie\CGF029';
dp = 'J:\AM152';

filelistcheck  = dir(fullfile(dp, '*.czi'));
filepaths      = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');

sliceinfo                = struct();
sliceinfo.filepaths      = filepaths;
sliceinfo.px_process     = 1.3;  % um
sliceinfo.px_register    = 20; % um
sliceinfo.slicethickness = 80; % um, slice thickness
sliceinfo                = getSliceInfo(sliceinfo);
%% we first generate the slice volume
slicevol = generateSliceVolume(sliceinfo);

%% (manual) reorder slices if needed
InteractiveSliceReorder(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (manual) flip slices if needed
SliceFlipper(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (auto) we align slices
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
sliceinfo.use_gpu  = true;
alignedvol = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

%% (auto) we find the atlas correspondence of our volume
% IN PROGRESS!!!!
% we load the atlas and get points on it

bulkAlignToAllen(sliceinfo)


%%
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
ichan      = find(contains(sliceinfo.channames, 'Cy3')); % We use the tdTomato channel!
sliceinfo.celldiam = 14; % expected cell diameter in um
celllocs = extractCellsFromSliceVolume(sliceinfo, ichan);