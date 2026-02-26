opts = struct();
%=========================================================================
% options to change
%--------------------------------------------------------------------------
% for naming
opts.mousename  = 'DK001'; 
% change for the folder that contains stitched tiff files
opts.datafolder      = 'D:\DATA\DK001'; 
opts.fproc      = fullfile('C:\DATA_sorted'); % where the processed volume is saved as a binary (fast SSD, at
% least 500 GB), will be deleted
% path to save results
opts.savepath   = fullfile(opts.datafolder, 'lightsuite'); 
%--------------------------------------------------------------------------
% some processing options
opts.tifftype           = 'channelperfile'; % can be planeperfile or channelperfile 
opts.pxsize             = [8.23 8.23 5]; % voxel size, xy and z, in um
opts.celldiam           = 14; % approximate cell size in um
opts.atlasres           = 10; % better keep this fixed for highest resolution, in um
opts.registres          = 20; % resolution to do the nonrigid registration, keep fixed, in um
opts.debug              = false; % toggle plotting (takes longer) for cell detections
opts.savecellimages     = false; % toggle saving of individual cell images
opts.thres_cell_detect  = [0.5 0.4]; % thresholds for detecting cells relative to background
opts.channelforcells    = 1; % channel to use for cell detection, leave empty ([]) for none
opts.channelforregister = 1; % channel to use for registration
opts.writetocsv         = false; % write cells to csv file
opts.augmentpoints      = false; % whether to augment user-defined points with automatic ones
opts.weight_usr_pts     = 0.2;   % set to zero for image-only information
%--------------------------------------------------------------------------
opts                   = readLightsheetOpts(opts);
%=========================================================================
%% (auto) main processing pipeline, preprocess and detect cells
opts = preprocessLightSheetVolume(opts);

%% (auto) initialize registration
opts = initializeRegistration(opts.savepath);

%% (manual) Refine control points
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
matchControlPoints_unified(opts);

%% (auto) perform registration
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
transform_params   = multiobjRegistration(opts, opts.weight_usr_pts, true);