
% folder which contains mouse subfolders
datafolderpath = 'D:\DATA';
mousename      = 'DK001';

dp        = fullfile(datafolderpath, sprintf('*%s*', mousename));
dp        = dir(dp);
dp        = fullfile(dp.folder, dp.name);
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

%% (manual) reorder, flip, center and discard slices if needed
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
opts.cpwt       = 0.4;
transformparams = registerSlicesToAtlas(opts);

%% (auto) apply registration to all color channels to generate registered volumes
transformparams = load(fullfile(sliceinfo.procpath, "transform_params.mat"));
sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo          = sliceinfo.sliceinfo;
generateRegisteredSliceVolume(sliceinfo, transformparams);

%% (auto) detect cells in slices
% Set ichan to a scalar or vector of channel indices to detect in.
% Any number of channels can be processed; results are saved per channel.
sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
ichans     = [2 3]; % e.g. [2] or [2 3] for multiple channels
sliceinfo.debug    = true;  % toggle plotting (takes longer) for detections
sliceinfo.celldiam = 14;    % expected cell diameter in um
sliceinfo.thresuse = single([0.75 0.4]); % thresholds in SBR (detection, expansion)
extractCellsFromSliceVolume(sliceinfo, ichans);

% you can use visualizeCellDetections to plot all the detections in sample
% space like this:
visualizeCellDetections(sliceinfo.procpath, Space = 'sample');

%% (auto) move cell detections in atlas space
% Results are saved as chan<NN>_cell_locations_atlas.mat in procpath.
transformparams = load(fullfile(sliceinfo.procpath, "transform_params.mat"));
for ci = 1:numel(ichans)
    curr_ichan = ichans(ci);
    celllocs   = load(fullfile(sliceinfo.procpath, ...
        sprintf('chan%02d_cell_locations_sample.mat', curr_ichan)));
    atlasptcoords = slicePointsToAtlas(celllocs.cell_locations, transformparams);
    fsavename = fullfile(sliceinfo.procpath, ...
        sprintf('chan%02d_cell_locations_atlas.mat', curr_ichan));
    save(fsavename, 'atlasptcoords')
end

% you can use visualizeCellDetections to plot all the detections in atlas
% space like this:
visualizeCellDetections(sliceinfo.procpath, Space = 'atlas');

%% (auto) write intensities per brain region to CSV
sliceinfo = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
sliceinfo = sliceinfo.sliceinfo;
generateSliceIntensitiesCSV(sliceinfo, 'writetocsv', true);

%% (auto) write cell counts per brain region to CSV
sliceCellCountsToCSV(sliceinfo.procpath, 'writetocsv', true);
