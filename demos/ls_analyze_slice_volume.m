
% folder which contains mouse subfolders
datafolderpath = 'D:\Histology\';
mousename      = 'AM130';

dp = fullfile(datafolderpath, sprintf('*%s*', mousename));
dp = dir(dp);
dp = fullfile(dp.folder, dp.name);
sliceinfo = parseSettingsFile(fullfile(dp, 'local_settings.txt'));

sliceinfo.mousename = mousename;

% --- Locate input files: CZI or a folder of TIF/TIFF images ---
filelistcheck = dir(fullfile(dp, '*.czi'));
if ~isempty(filelistcheck)
    % CZI mode: each file may contain multiple scenes (slices)
    filepaths = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
else
    % TIF mode: one TIF file per slice, read directly (no scene selection).
    %   Multi-channel slices should be stored as multi-page TIFs (one page
    %   per channel). Single-channel TIFs are also supported.
    filelistcheck = dir(fullfile(dp, '*.tif'));
    if isempty(filelistcheck)
        filelistcheck = dir(fullfile(dp, '*.tiff'));
    end
    filepaths = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');

    % Pixel size and channel names cannot be read from TIF metadata and
    % must be specified manually before calling getSliceInfo.
    sliceinfo.pxsize    = [0.65 0.65];           % um/pixel — set to your value
    sliceinfo.channames = {'DAPI', 'Cy3', 'Cy5'}; % set to your channel names
end

sliceinfo.filepaths = filepaths;
sliceinfo           = getSliceInfo(sliceinfo);

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

% we can plot the cells in atlas space
nrand = min(size(atlasptcoords,1), 1e5);
iplot = randperm(size(atlasptcoords,1),nrand);
close all;
plotBrainGrid; hold on;
scatter3(atlasptcoords(iplot,2),atlasptcoords(iplot,3),atlasptcoords(iplot,1),2,'filled','MarkerFaceAlpha',0.5)
