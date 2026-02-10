
mouselist = {'YX020','YX021','YX026','YX027','YX028', 'DK031', 'DK025','DK026','DK027','DK028'};
for ibrain = 1:numel(mouselist)

    opts = struct();
    %--------------------------------------------------------------------------
    % this is the stitched file path
    opts.mousename  = mouselist{ibrain};
    disp(opts.mousename)
    dpfind          = fullfile('F:\imaging', sprintf('%s*',opts.mousename), '**','*.tif');
    allfilesfind    = dir(dpfind);
    if isempty(allfilesfind)
        continue
    end
    [allun, ~, iun] = unique({allfilesfind(:).folder}');
    [~, ifolder]    = max(accumarray(iun, 1, [numel(allun) 1], @sum));
    opts.datafolder = allun{ifolder};
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
    %% (auto) Preprocess high-res volume and extract cells
    [~, opts]  = preprocessColmVolumeNew(opts);
    opts.debug = false; % toggle plotting (takes longer) for detections
    opts.savecellimages = true;
    gpuDevice(1);
    cell_locations = extractCellsFromVolumeNew(opts);
    delete(opts.fproc);
    %%
    transform_params = load(fullfile(opts.savepath, 'transform_params.mat'));
    
    % load and transform background volume
    backvolume                  = readDownStack(fullfile(opts.savepath, 'sample_register_20um.tif'));
    [backvolareas, areainds]    = backVolumeToAtlas(backvolume, transform_params);
    fsavename                   = fullfile(opts.savepath, 'background_volume_areas.mat');
    save(fsavename, 'backvolareas', 'areainds') 

    % load and transform cell locations
    celllocs         = load(fullfile(opts.savepath,'cell_locations_sample.mat'), 'cell_locations');
    atlasptcoords    = volumePointsToAtlas(celllocs.cell_locations, transform_params);
    fsavename        = fullfile(opts.savepath, 'cell_locations_atlas.mat');
    save(fsavename, 'atlasptcoords') 
end