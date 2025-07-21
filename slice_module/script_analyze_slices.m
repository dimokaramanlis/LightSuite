
% folder which contains mouse subfolders
datafolderpath = 'J:\'; %'D:\example_charlie'; %'J:\'; % 
mousename      = 'AM130';%'AM147';%'CGF028'; %'AM130'; 'CGF027';

dp = fullfile(datafolderpath, sprintf('*%s*', mousename));
dp = dir(dp);
dp = fullfile(dp.folder, dp.name);
sliceinfo = parseSettingsFile(fullfile(dp, 'local_settings.txt'));

sliceinfo.mousename   = mousename;
filelistcheck         = dir(fullfile(dp, '*.czi'));
filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
sliceinfo.filepaths   = filepaths;
sliceinfo             = getSliceInfo(sliceinfo);
%% (auto) we first generate the slice volume
slicevol = generateSliceVolume(sliceinfo, sliceinfo.regchan);

%% (manual) reorder, flip and discard slices if needed
SliceOrderEditor(sliceinfo.volorder)
generateReordedVolume(sliceinfo);

%% (auto) we align slices and initialize registration
sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo          = sliceinfo.sliceinfo;
sliceinfo          = copyStructBtoA(sliceinfo, settings);
alignedvol         = alignSliceVolume(sliceinfo.slicevol, sliceinfo);
%% (manual) determine cutting angle gui if you are not happy with the original estimation
opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
determineCuttingAngleGUI(opts)
%% (manual) match control points to determine cutting angle and gaps

% !!! The control point selection is currently tied to the initial
% registration. Don't start before checking the diagnostic plots and the
% inspection volume!!!

opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
matchControlPointsInSlices(opts)
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
cellsout         = cellsout.cell_locations;
transformparams  = load(fullfile(sliceinfo.procpath, "transform_params.mat"));
atlasptcoords    = slicePointsToAtlas(cellsout, transformparams);
fsavename        = fullfile(sliceinfo.procpath, 'cell_locations_atlas.mat');
save(fsavename, 'atlasptcoords') 

% nrand = min(size(atlasptcoords,1), 1e5);
% iplot = randperm(size(atlasptcoords,1),nrand);
% close all;
% plotBrainGrid; hold on;
% scatter3(atlasptcoords(iplot,2),atlasptcoords(iplot,3),atlasptcoords(iplot,1),2,'filled','MarkerFaceAlpha',0.5)