
% folder which contains mouse subfolders
datafolderpath = 'J:\'; % 'D:\example_charlie'; %
sliceinfo                = struct();
sliceinfo.mousename      = 'AM130';%'AM147';%'CGF028'; %'AM130'; 

%slice thickness in um, 150 for thicker slices
sliceinfo.slicethickness =  80; 
%size for processing slices in um, 1.25 for cell detection, 5 for other signals
sliceinfo.px_process     = 1.25; 

dp = fullfile(datafolderpath, sprintf('*%s*', sliceinfo.mousename));
dp = dir(dp);
dp = fullfile(dp.folder, dp.name);

filelistcheck         = dir(fullfile(dp, '*.czi'));
filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
sliceinfo.filepaths   = filepaths;
sliceinfo.px_register = 20; % um
sliceinfo.px_atlas    = 10; % um
sliceinfo.atlasaplims = [180 1079]; % these limits determine where we look in the atlas
sliceinfo             = getSliceInfo(sliceinfo);
%% (auto) we first generate the slice volume
sliceinfo.denoisedapi   = true; % activate if your dapi looks like crap
slicevol = generateSliceVolume(sliceinfo);

%% (manual) reorder, flip and discard slices if needed
SliceOrderEditor(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (auto) we align slices and initialize registration
sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo          = sliceinfo.sliceinfo;
sliceinfo.use_gpu  = true;
sliceinfo.regmedianfilt = true; % activate if you want to register with cell channel 
alignedvol         = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

%% (manual) match control points to determine cutting angle and gaps

% !!! The control point selection is currently tied to the initial
% registration. Don't start before checking the diagnostic plots and the
% inspection volume!!!

opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
matchControlPointsInSliceVolume(opts)


%% (auto) refine registation with control points and elastix
opts            = load(fullfile(sliceinfo.procpath, "regopts.mat"));
transformparams = registerSlicesToAtlas(opts);

%% (auto) apply registration to all color channels to generate registered volumes
transformparams = load(fullfile(sliceinfo.procpath, "transform_params.mat"));
sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo          = sliceinfo.sliceinfo;
generateRegisteredSliceVolume(sliceinfo, transformparams);

%% (auto) detect cells in slices
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
ichan      = find(contains(sliceinfo.channames, 'Cy3')); % We use the tdTomato channel!
sliceinfo.debug    = true; % toggle plotting (takes longer) for detections
sliceinfo.celldiam = 14; % expected cell diameter in um
sliceinfo.thresuse = single([0.75 0.4]); % thresholds in SBR units (first for detection and then for expansion)
celllocs = extractCellsFromSliceVolume(sliceinfo, ichan);

%% (auto) move cell detections in atlas space
cellsout         = load(fullfile(sliceinfo.procpath,'cell_locations_sample.mat'));
transformparams  = load(fullfile(sliceinfo.procpath, "transform_params.mat"));
atlasptcoords    = slicePointsToAtlas(cellsout, transformparams);
fsavename        = fullfile(opts.savepath, 'cell_locations_atlas.mat');
save(fsavename, 'atlasptcoords') 