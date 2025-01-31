opts = struct();
%--------------------------------------------------------------------------
% this is the stitched file path
opts.mousename  = 'DK027';
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
opts.celldiam   = 10; % in um
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
peakvalsextract = preprocessColmVolumeBatches(opts);
% irand = randperm(size(peakvalsextract,1), 100000);
% scatter3(peakvalsextract(irand,1)*opts.pxsize(1),...
%     peakvalsextract(irand,2)*opts.pxsize(2),...
%     peakvalsextract(irand,3)*opts.pxsize(3),1, peakvalsextract(irand,5))
% axis equal; ax = gca; ax.ZDir = 'reverse';
opts           = initializeRegistration(opts);
%%
% we perform initial registration to ease control point selection
opts.regvolpath = fullfile('C:\DATA_sorted', sprintf('%s_temporary_reg_volume.dat', opts.mousename));
opts.regvolsize = ceil([opts.Ny opts.Nx opts.Nz].*opts.pxsize/opts.registres);

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
transform_params = multiobjRegistration(opts, 1);


%%
celllocs      = load(fullfile(opts.savepath,'cell_locations_sample.mat'));
atlasptcoords = volumePointsToAtlas(celllocs.cell_locations, transform_params);
fsavename     = fullfile(opts.savepath, 'cell_locations_atlas.mat');
save(fsavename, 'atlasptcoords')
%%

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));

st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
%%
cleanatlaspts                    = sanitizeCellCoords(atlasptcoords, av);
[areacounts, areavols, idcareas] = groupCellsIntoLeafRegions(cleanatlaspts, av);
stmat                            = atlasTreeToMat(st);
Nmax = numel(areacounts);
%%
groupstrs       = lower(table2cell(st(1:Nmax,4)));
areavolsuse     = areavols; % in mm3
cells_per_group = areacounts;
irem = cells_per_group<2;
cells_per_group(irem) = [];
areavolsuse(irem) = [];
stgroups = stmat(1:Nmax, :);

% % move one up
% lastleaf   = findfirst(~isnan(stgroups),2,1,'last');
% prevbranch = lastleaf - 1;
% idxprev    = sub2ind(size(stgroups), find(prevbranch>0), prevbranch(prevbranch>0));
% 
% 
% idxprev - 1

groupstrs(irem) = [];




figure;
Nshow = 50;
cell_densities = cells_per_group./areavolsuse;
% cell_densities = cells_per_group;

[~, isortc] = sort(cell_densities, 'descend');
bar(groupstrs(isortc(1:Nshow)), cell_densities(isortc(1:Nshow)))

ivis = contains(groupstrs(isortc),'primary motor')|...
    contains(groupstrs(isortc),'primary auditory')|...
    contains(groupstrs(isortc),'motor')|...
    contains(groupstrs(isortc),'primary visual');


ivis = contains(groupstrs(isortc),'primary motor')|...
    contains(groupstrs(isortc),'secondary motor')|...
    contains(groupstrs(isortc),'prelimbic')|...
    contains(groupstrs(isortc),'infralimbic')|...
    contains(groupstrs(isortc),'cingulate');

bar(groupstrs(isortc(ivis)), cell_densities(isortc(ivis)))
ylabel('TRAP cell density (cells/mm^3)')
%%

avplot = ndSparse.build(cleanatlaspts(:,[2 1 3]), 1,size(av));
avplot = single(full(avplot));
avgauss = imgaussfilt3(avplot, 10);

%%
% extract 3d maxima over batches to generate candidate list
% running over batches again
% for each batch (+ overlap a bit), get brightest point, fit 3D Gaussian,
% remove point from list and all points that are within 2sigma of that
% Gaussian, and remove Gaussian from volume

%%
% for background imadjust(imflatfield(backsignal, 1000))

