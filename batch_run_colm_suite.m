 opts = struct();
%--------------------------------------------------------------------------
% this is the stitched file path
opts.mousename  = 'YX017';
dpfind          = fullfile('F:\imaging', sprintf('%s*',opts.mousename), '**','*.tif');
allfilesfind    = dir(dpfind);
[allun, ~, iun] = unique({allfilesfind(:).folder}');
[~, ifolder]    = max(accumarray(iun, 1, [numel(allun) 1], @sum));
opts.datafolder = allun{ifolder};
% opts.datafolder = 'F:\imaging\DK025_uncompressed\VW0\Stitched\RES(10340x7449x1889)\000000\000000_0-1550';
% opts.datafolder = 'S:\ElboustaniLab\#SHARE\Data\DK031\Anatomy\colm\561_RES(11672x8464x1581)';
%--------------------------------------------------------------------------
% this is where the processed volume is saved as a binary (fast SSD, at least 500 GB)
opts.fproc      = fullfile('C:\DATA_sorted', sprintf('%s_fullbrain_dff.dat', opts.mousename));
opts.savepath   = fullfile('D:\DATA_folder\Mice', opts.mousename, 'Anatomy');
opts.pxsize     = [1.44 1.44 5]; % in um
opts.celldiam   = 14; % in um
opts.atlasres   = 10; % in um
opts.registres  = 20; % in um
% some processing options
opts.maxdff     = 12;
opts            = readLightsheetOpts(opts);
%%
% start with preprocessing and extract candidate cells
% this part is slow
[backvol, opts] = preprocessColmVolume(opts);
% let's make sure the GPU can support our endeavor
gpuDevice(1);
peakvalsextract = extractCellsFromVolume(opts.fproc, opts);
delete(opts.fproc)
irand = randperm(size(peakvalsextract,1), 80000);
scatter3(peakvalsextract(irand,1)*opts.pxsize(1),...
    peakvalsextract(irand,2)*opts.pxsize(2),...
    peakvalsextract(irand,3)*opts.pxsize(3),1)
axis equal; ax = gca; ax.ZDir = 'reverse';
opts           = initializeRegistration(opts);
%%
% SELECT CONTROL POINTS
% In the function, the sample volume goes through a rigid transform to
% ease control point selection
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
matchControlPoints(opts);
%%
% MAIN REGISTRATION FUNCTION
% we have to move selected control points back in sample space
usemultistep     = true;
transform_params = multiobjRegistration(opts, 0.2, usemultistep);

%%
transform_params = load(fullfile(opts.savepath, 'transform_params.mat'));

% load and transform background volume
backvolume                  = readDownStack(fullfile(opts.savepath, 'sample_register_20um.tif'));
[backvolareas, areainds]    = backVolumeToAtlas(backvolume, transform_params);
fsavename                   = fullfile(opts.savepath, 'background_volume_areas.mat');
save(fsavename, 'backvolareas', 'areainds') 
plot(backvolareas(:,1), backvolareas(:,2),'.',[0 20],[0 20])
axis equal; xlim([0 20]); ylim([0 20])
title(opts.mousename)

%%
% load and transform cell locations
celllocs         = load(fullfile(opts.savepath,'cell_locations_sample.mat'));
atlasptcoords    = volumePointsToAtlas(celllocs.cell_locations, transform_params);
fsavename        = fullfile(opts.savepath, 'cell_locations_atlas.mat');
save(fsavename, 'atlasptcoords') 

% moveToAtlasTest(celllocs.cell_locations,backvolume, transform_params);
%%

% celllocs         = load(fullfile(opts.savepath,'cell_locations_atlas.mat'));
% atlasptcoords = celllocs.atlasptcoords;

nrand = min(size(atlasptcoords,1), 1e5);
iplot = randperm(size(atlasptcoords,1),nrand);
plotBrainGrid; hold on;
scatter3(atlasptcoords(iplot,2),atlasptcoords(iplot,3),atlasptcoords(iplot,1),2,'filled','MarkerFaceAlpha',0.5)
