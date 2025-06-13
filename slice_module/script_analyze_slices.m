
datafolderpath = 'D:\example_charlie';

sliceinfo                = struct();
sliceinfo.mousename      = 'CGF028';
sliceinfo.slicethickness = 150; % um, slice thickness
sliceinfo.px_process     = 5;  % um

dp = fullfile(datafolderpath, sprintf('*%s*', sliceinfo.mousename));
dp = dir(dp);
dp = fullfile(dp.folder, dp.name);

filelistcheck         = dir(fullfile(dp, '*.czi'));
filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
sliceinfo.filepaths   = filepaths;
sliceinfo.px_register = 20; % um
sliceinfo.px_atlas    = 10; % um
sliceinfo.atlasaplims = [180 1079];
sliceinfo             = getSliceInfo(sliceinfo);
%% (auto) we first generate the slice volume
slicevol = generateSliceVolume(sliceinfo);

%% (manual) reorder slices if needed
InteractiveSliceReorder(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (manual) flip slices if needed
SliceFlipper(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (auto) we align slices and initialize registration
sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo          = sliceinfo.sliceinfo;
sliceinfo.use_gpu  = true;
alignedvol         = alignSliceVolume(sliceinfo.slicevol, sliceinfo);
%% (manual) match control points
opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
matchControlPointsInSlices(opts)

%%
optsfin = registerSlicesToAtlas(opts);

%% (auto) we find the atlas correspondence of our volume
% IN PROGRESS!!!!
% we load the atlas and get points on it
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
sliceinfo.use_gpu  = true;

bulkAlignToAllen(sliceinfo)
%%
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
writeVolumeForDeepSlice(sliceinfo, 'DAPI')

%%
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
ichan      = find(contains(sliceinfo.channames, 'Cy3')); % We use the tdTomato channel!
sliceinfo.celldiam = 14; % expected cell diameter in um
celllocs = extractCellsFromSliceVolume(sliceinfo, ichan);