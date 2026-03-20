function varargout = transformPointsToAtlas(input_data, varargin)
%TRANSFORMPOINTSTOATLAS Transforms cell coordinates to atlas space and counts them.
%
%   FINALPTS = TRANSFORMPOINTSTOATLAS(INPUT_PTS, 'transform_params', TRSTRUCT) 
%   transforms the N x M array of INPUT_PTS directly. N is the number of 
%   detected points, and M must be >= 3. The first 3 columns represent the 
%   [x, y, z] spatial coordinates in the original sample space. Any additional 
%   columns (e.g., intensity, equivalent diameter, ellipticity) are purely 
%   descriptive and will be carried over to the output unmodified. TRSTRUCT 
%   must contain the transformation parameters.
%
%   TRANSFORMPOINTSTOATLAS(SAVEPATH) searches the directory SAVEPATH for 
%   '*cell_locations_sample.mat' files, transforms the points, saves them 
%   as '*_atlas.mat' inside the 'volume_registered' directory, and calculates 
%   cell counts per brain region using the Allen Brain Atlas.
%
%   TRANSFORMPOINTSTOATLAS(..., 'writetocsv', VAL) specifies whether to 
%   write the transformed statistics (and point locations) to CSV files.
%
%   Inputs:
%       input_data           - (char/string) Directory path OR 
%                              (numeric) N x M array of points (N points, 
%                              M dimensions where cols 1:3 are x, y, z 
%                              in sample space, and cols 4:end are single 
%                              number descriptors like intensity/diameter).
%
%   Optional Name-Value Parameters:
%       transform_params     - (struct) Transformation structure. Required 
%                              if input_data is a numeric array.
%       writetocsv           - (logical) If true, writes data to CSV. 
%
%   Outputs:
%       finalpts             - (single/double) Array of transformed points. 
%                              Returned if an output variable is requested.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1. Parse Inputs and establish defaults
%--------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'input_data', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'transform_params', [], @isstruct);
addParameter(p, 'writetocsv', [], @(x) islogical(x) || isscalar(x));

parse(p, input_data, varargin{:});

if ischar(input_data) || isstring(input_data)
    opMode   = 'folder';
    savepath = char(input_data);
else
    opMode   = 'points';
    inputpts = input_data;
end

%--------------------------------------------------------------------------
% 2. Execute based on Operation Mode
%--------------------------------------------------------------------------
if strcmp(opMode, 'folder')
    fprintf('Running in folder mode. Target: %s\n', savepath);
    
    % Load configuration and transforms
    trfile  = fullfile(savepath, 'transform_params.mat');
    optfile = fullfile(savepath, 'regopts.mat');
    
    if ~exist(trfile, 'file') || ~exist(optfile, 'file')
        error('transform_params.mat or regopts.mat not found in the folder.');
    end
    
    trstruct   = load(trfile);
    loadedOpts = load(optfile);
    opts       = loadedOpts.opts;
    
    % Determine writetocsv fallback
    writetocsv = false;
    if ~isempty(p.Results.writetocsv)
        writetocsv = p.Results.writetocsv;
    elseif isfield(opts, 'writetocsv')
        writetocsv = opts.writetocsv;
    end
    
    % Find point files
    files = dir(fullfile(savepath, '*cell_locations_sample.mat'));
    if isempty(files)
        fprintf('No valid cell location files found to process.\n');
        if nargout > 0; varargout{1} = []; end
        return;
    end
    
    % Setup registration output directory
    registerpath = fullfile(savepath, 'volume_registered');
    if ~exist(registerpath, 'dir')
        mkdir(registerpath);
    end
    
    % Load Allen Atlas and Parcellation info once
    fprintf('Loading Allen Atlas and parcellation data...\n');
    atlasData = loadAtlasData();
    groupinds = atlasData.areaidx;
    
    % Process each file
    all_final_pts = cell(length(files), 1);
    for i = 1:length(files)
        currfile = fullfile(files(i).folder, files(i).name);
        fprintf('Processing %s... \n', files(i).name); savetic = tic;
        
        dat = load(currfile, 'cell_locations');
        if ~isfield(dat, 'cell_locations')
            fprintf('  -> Skipping: ''cell_locations'' not found.\n');
            continue;
        end
        
        % 1. Transform Points
        atlasptcoords = coreTransform(dat.cell_locations, trstruct);
        
        % 2. Sanitize and remove invalid points
        [cleanatlaspts, badpts] = sanitizeCellCoords(atlasptcoords, atlasData.av);
        
        % 3. Group into brain regions and calculate stats
        [areacounts, areavols, catids] = groupCellsIntoLeafRegions(...
            cleanatlaspts, atlasData.av, groupinds);
        
        % 4. Save locations in the volume_registered folder
        outname   = strrep(files(i).name, '_sample.mat', '_atlas.mat');
        fsavename = fullfile(registerpath, outname);
        save(fsavename, 'atlasptcoords', 'cleanatlaspts', 'badpts', 'catids');

        if writetocsv
            loc_csvname = strrep(outname, '.mat', '.csv');
            writematrix(atlasptcoords, fullfile(registerpath, loc_csvname));
        end
        
        % 5. Extract channel number for stats naming
        tokens = regexp(files(i).name, 'chan_(\d+)_', 'tokens');
        if isempty(tokens)
            ichan = i; % Fallback if regex fails
        else
            ichan = str2double(tokens{1}{1});
        end
        
        % 6. Save statistics
        saveCellStats(registerpath, ichan, areacounts, areavols, atlasData, writetocsv);
        
        all_final_pts{i} = atlasptcoords;
        fprintf('  -> Done! Channel %d completed in %2.2f s.\n', ichan, toc(savetic));
    end
    
    if nargout > 0
        if length(all_final_pts) == 1
            varargout{1} = all_final_pts{1};
        else
            varargout{1} = all_final_pts;
        end
    end
    fprintf('All points successfully registered and tabulated!\n');
    
else
    % Points mode (direct transformation without statistics)
    trstruct = p.Results.transform_params;
    if isempty(trstruct)
        error('''transform_params'' structure must be provided when passing points directly.');
    end
    
    finalpts = coreTransform(inputpts, trstruct);
    if nargout > 0
        varargout{1} = finalpts;
    end
end

end

%==========================================================================
% Local Helper Functions
%==========================================================================

function finalpts = coreTransform(inputpts, trstruct)
    regsize_mm = trstruct.atlasres * 2 * 1e-3;
    
    pts = (inputpts(:, 1:3) - 1) .* trstruct.ori_pxsize * 1e-3; 
    pts = pts(:, [2 1 3]);
    
    if isfield(trstruct, 'ori_size')
        phys_size_yxz = (trstruct.ori_size - 1) .* trstruct.ori_pxsize([2 1 3]) * 1e-3;
    else
        error('trstruct must contain ''ori_size'' for axis flips.');
    end
    
    perm_order   = abs(trstruct.how_to_perm);
    pts_permuted = zeros(size(pts), 'like', pts);
    
    for dim = 1:3
        orig_dim = perm_order(dim);
        pts_permuted(:, dim) = pts(:, orig_dim);
        if trstruct.how_to_perm(dim) < 0
            dim_max_size = phys_size_yxz(orig_dim);
            pts_permuted(:, dim) = dim_max_size - pts_permuted(:, dim);
        end
    end
    
    pts = pts_permuted(:, [2 1 3]) / regsize_mm;
    
    Dfield = transformix([], trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    Dfield = permute(Dfield, [2 3 4 1]) / regsize_mm; 
    [Sx, Sy, Sz, ~] = size(Dfield);
    
    Xgv = 1:Sx; Ygv = 1:Sy; Zgv = 1:Sz;
    dx = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,1), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dy = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,2), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dz = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,3), pts(:,1), pts(:,2), pts(:,3), 'linear');
    
    interpolated_displacements = -[dx, dy, dz];
    interpolated_displacements(isnan(interpolated_displacements)) = 0;
    
    finalpts = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward(pts + interpolated_displacements);
    finalpts = [finalpts, inputpts(:, 4:end)];
end

function atlasData = loadAtlasData()
    % Loads the annotation volume and parcellation CSV info
    allen_atlas_path        = fileparts(which('annotation_10.nii.gz'));
    av                      = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
    allen_atlas_parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));

    parcelinfo       = readtable(fullfile(allen_atlas_parcel_path, 'parcellation_to_parcellation_term_membership.csv'));
    
    substridx        = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib]    = unique(parcelinfo.parcellation_index(substridx));
    namessub         = parcelinfo.parcellation_term_name(substridx);
    
    stridx           = strcmp(parcelinfo.parcellation_term_set_name, 'structure');
    [~, ibstr]       = unique(parcelinfo.parcellation_index(stridx));
    namesstruct      = parcelinfo.parcellation_term_name(stridx);
    
    dividx           = strcmp(parcelinfo.parcellation_term_set_name, 'division');
    [~, ibdiv]       = unique(parcelinfo.parcellation_index(dividx));
    namesdiv         = parcelinfo.parcellation_term_name(dividx);

    % Pack into a clean struct
    atlasData.av          = av;
    atlasData.parcelinfo  = parcelinfo;
    atlasData.areaidx     = areaidx;
    atlasData.namessub    = namessub(ib);
    atlasData.namesstruct = namesstruct(ibstr);
    atlasData.namesdiv    = namesdiv(ibdiv);
end

function saveCellStats(registerpath, ichan, areacounts, areavols, atlasData, writetocsv)
    % Saves the cell counts mirroring generateRegisteredBrainVolumes
    fmatname = fullfile(registerpath, sprintf('chan%02d_cellcounts.mat', ichan));
    areaidx  = atlasData.areaidx;
    save(fmatname, 'areacounts', 'areaidx', 'areavols');

    if writetocsv
        currtable = array2table([areaidx, areacounts, areavols], ...
            'VariableNames', ...
            {'parcellation_index', 'RightSideCount', 'LeftSideCount', 'TotalVolume[mm3]'});
            
        currtable = addvars(currtable, atlasData.namessub, atlasData.namesstruct, atlasData.namesdiv, ...
            'NewVariableNames', {'name', 'structure', 'division'}, 'Before', 'parcellation_index');
            
        fsavename = fullfile(registerpath, sprintf('chan%02d_cellcounts.csv', ichan));
        writetable(currtable, fsavename);
    end
end