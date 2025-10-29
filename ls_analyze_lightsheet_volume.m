opts = struct();
%=========================================================================
% options to change
%--------------------------------------------------------------------------
% for naming
opts.mousename  = 'DK001'; 
% change for the folder that contains stitched tiff files
tifffolder      = 'C:\DATA\DK001'; 
% this is where the processed volume is saved as a binary (fast SSD, at
% least 500 GB), will be deleted
opts.fproc      = fullfile('C:\DATA_sorted', sprintf('%s_fullbrain_dff.dat', opts.mousename));
% path to save results
opts.savepath   = fullfile(tifffolder, 'lightsuite'); 
%--------------------------------------------------------------------------
% some processing options
opts.maxdff         = 12; % parameter for preprocessing
opts.pxsize         = [1.44 1.44 5]; % voxel size, xy and z, in um
opts.celldiam       = 14; % approximate cell size in um
opts.atlasres       = 10; % better keep this fixed for highest resolution, in um
opts.registres      = 20; % resolution to do the nonrigid registration, keep fixed, in um
opts.debug          = false; % toggle plotting (takes longer) for cell detections
opts.savecellimages = false; % toggle saving of individual cell images
%--------------------------------------------------------------------------
dpfind          = fullfile(tifffolder, '**','*.tif');
allfilesfind    = dir(dpfind);
[allun, ~, iun] = unique({allfilesfind(:).folder}');
[~, ifolder]    = max(accumarray(iun, 1, [numel(allun) 1], @sum));
opts.datafolder = allun{ifolder};
%--------------------------------------------------------------------------
opts            = readLightsheetOpts(opts);
%=========================================================================
%% main processing pipeline, preprocess and detect cells
[backvol, opts] = preprocessColmVolume(opts);
% let's make sure the GPU can support our endeavor
gpuDevice(1);
peakvalsextract = extractCellsFromVolume(opts.fproc, opts);
delete(opts.fproc)
%% test cell detections
irand = randperm(size(peakvalsextract,1), 80000);
scatter3(peakvalsextract(irand,1)*opts.pxsize(1),...
    peakvalsextract(irand,2)*opts.pxsize(2),...
    peakvalsextract(irand,3)*opts.pxsize(3),1)
axis equal; ax = gca; ax.ZDir = 'reverse';
opts           = initializeRegistration(opts);
%% manually curate registration
% In the function, the sample volume goes through a rigid transform to
% ease control point selection
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
matchControlPoints(opts);
%% perfom full registration (sanity check)
% we have to move selected control points back in sample space
usemultistep     = true;
transform_params = multiobjRegistration(opts, 0.2, usemultistep);
%% load and transform background volume

transform_params = load(fullfile(opts.savepath, 'transform_params.mat'));
backvolume                  = readDownStack(fullfile(opts.savepath, 'sample_register_20um.tif'));
[backvolareas, areainds]    = backVolumeToAtlas(backvolume, transform_params);
fsavename                   = fullfile(opts.savepath, 'background_volume_areas.mat');
save(fsavename, 'backvolareas', 'areainds') 
%% load and transform cell locations
celllocs         = load(fullfile(opts.savepath,'cell_locations_sample.mat'));
atlasptcoords    = volumePointsToAtlas(celllocs.cell_locations, transform_params);
fsavename        = fullfile(opts.savepath, 'cell_locations_atlas.mat');
save(fsavename, 'atlasptcoords') 
%% plot detected cell locations in atlas space (sanity check)

celllocs         = load(fullfile(opts.savepath,'cell_locations_atlas.mat'));
atlasptcoords = celllocs.atlasptcoords;

nrand = min(size(atlasptcoords,1), 1e5);
iplot = randperm(size(atlasptcoords,1),nrand);
plotBrainGrid; hold on;
scatter3(atlasptcoords(iplot,2),atlasptcoords(iplot,3),atlasptcoords(iplot,1),2,...
    'filled','MarkerFaceAlpha',0.5)
