opts = struct();
%=========================================================================
% options to change
%--------------------------------------------------------------------------
% for naming
opts.mousename  = 'test';
% change for the folder that contains the stitched tiff files
opts.datafolder = 'F:\DATAcord';
opts.fproc      = fullfile('C:\DATA_sorted'); % where the processed volume is saved as a binary (fast SSD),
% will be deleted
% path to save results
opts.savepath   = fullfile(opts.datafolder, 'lightsuite');
%--------------------------------------------------------------------------
% some processing options
opts.tifftype           = 'channelperfile'; % can be planeperfile or channelperfile
opts.pxsize             =  [2.75 2.75 5]; % voxel size, xy and z, in um
opts.registres          = 20; % resolution to do the nonrigid registration, keep fixed, in um
opts.usegpu             = true; % activate if you have a GPU, used for cell detection
% cell detection parameters
opts.debug              = true; % toggle plotting (takes longer) for cell detections
opts.savecellimages     = true; % toggle saving of individual cell images
opts.celldiam           = 10; % approximate cell size in um
opts.thres_cell_detect  = [0.5 0.3]; % thresholds for detecting cells relative to background, first should be bigger than second
opts.channelforcells    = []; % channel to use for cell detection, leave empty ([]) for none
opts.writetocsv         = true; % write results to csv files
%  registration
opts.channelforregister = 2; % channel to use for registration
%--------------------------------------------------------------------------
opts                    = readLightsheetOpts(opts);
%=========================================================================
%% (auto) main processing pipeline, preprocess and detect cells
% same function as for a brain: it downsamples every channel to the
% registration resolution and, if opts.channelforcells is set, detects cells
opts = preprocessLightSheetVolume(opts);
%%
% you can also use visualizeCellDetections to plot all the detections in sample
% space like this:
visualizeCellDetections(opts.savepath, Space = 'sample');
%% (auto) initialize registration, part 1: orient, segment, crop, predict
% a cord is a long twisted tube cut out of a larger block, so before the atlas
% can be fitted the pipeline works out which axis is rostrocaudal, which end is
% rostral, where the brain starts, and - for every remaining slice - where the
% section sits and how it is rotated. Check cord_automatic_centerline.png in the
% save folder before moving on.
opts = prepareCordSampleForRegistration(opts.savepath);

%% (manual) correct the predicted centre line where it is wrong
% the GUI fits one spline through your clicks and the prediction above, dropping
% the prediction within opts.userexclusionradius slices of anything you click -
% so a click owns its stretch of cord outright and only the slices where the
% automatic answer is off need attention. The dashed prediction in the side
% plots disappears where you have taken over. If the whole cord comes out
% rotated by 180 degrees, press 'f' once. 'l' sets the regularisation and the
% reach of a click. Press 's' to save.
opts = loadRegOpts(opts.savepath);
spinal_cord_aligner(opts);

%% (auto) initialize registration, part 2: straighten and fit the atlas
% straightens the cord slice by slice, then stretches and affinely fits the
% atlas onto it. Check registration_initial_affine.png.
opts = initializeCordRegistration(opts.savepath);

%% (manual) add control points
opts = loadRegOpts(opts.savepath);
matchControlPoints_unified(opts);

%% (auto) perform nonlinear registration (b-spline)
opts = loadRegOpts(opts.savepath);
% options for registration
opts.weight_usr_pts        = 0.5;  % weight of user-defined points, set to zero for image-only information
opts.bspline_spatial_scale = 0.96; % in mm, how much you allow the bspline to bend
transform_params = multiobjCordRegistration(opts.savepath, opts.weight_usr_pts, ...
    'bspline_spatial_scale', opts.bspline_spatial_scale);

%% (auto) apply registration to all volume channels
% returns the volumes in atlas space, on the native atlas grid, so the
% annotation from loadSpinalCordAtlas indexes them voxel for voxel. Each region
% is also summarized per spinal segment by its median voxel intensity; pass
% 'areafun' to use something else (@mean, @std, or any handle taking a vector
% and returning one number). Every region is named the way the brain outputs
% name theirs - substructure, structure and division - so the CSVs can be read
% without looking ids up in the atlas.
[atlasvol, areastats] = generateRegisteredCordVolume(opts.savepath, ...
    'writetocsv', true, 'saveregisteredvolume', true);

%% (auto) apply registration to cell detections
% same file formats as for a brain: LightSuite .mat detections, the same array
% as .csv, and ImageJ "Cell Counter" .xml marker files. Counts come out per
% atlas region and per spinal segment, named at all three levels. Aggregate them
% with reorganizeSpinalCordAreas(counts, [], volumes, parcelinfo, areaidx,
% 'structure') - or 'division', or 'substructure'.
transformCordPointsToAtlas(opts.savepath, 'writetocsv', true);

% Points counted outside LightSuite go through the same transform, as long as
% they name their channel:
%   transformCordPointsToAtlas('D:\data\cells.xml', 'savepath', opts.savepath, ...
%       'channel', [2 3], 'writetocsv', true);   % Type 1 -> chan 2, Type 2 -> chan 3
%
% To drop detection artifacts with the CNN classifier, pass a trained network
% (see demos\trainClassificationNetwork.m). This requires opts.savecellimages =
% true during detection:
%   transformCordPointsToAtlas(opts.savepath, 'writetocsv', true, ...
%       'network', 'D:\nets\20260602_CellClassifierNet.mat');
