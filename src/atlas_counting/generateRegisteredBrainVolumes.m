function allmedians = generateRegisteredBrainVolumes(savepath, varargin)
%GENERATEREGISTEREDBRAINVOLUMES Registers brain volumes and calculates intensities.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(SAVEPATH) loads transformation 
%   parameters and options from SAVEPATH, applies the registration to the 
%   volumes, and calculates the background fluorescence in atlas coordinates.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(..., 'writetocsv', VAL) 
%   specifies whether to write the parcellation intensities to a CSV file.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(..., 'saveregisteredvolume', VAL)
%   specifies whether to save the registered volumes to disk.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(..., 'areafun', FUN) summarizes
%   the voxels of each area with FUN instead of the median. FUN must be a
%   function handle taking a column vector and returning a scalar - @median
%   (default), @mean and @std are the obvious choices, but anything of that
%   shape works (e.g. @(x) quantile(x, 0.9)). The same function is used for the
%   out-of-brain background level, so the two stay comparable.
%
%   Inputs:
%       savepath             - (char/string) Directory containing the volume
%                              files, 'transform_params.mat', and 'regopts.mat'.
%
%   Optional Name-Value Parameters:
%       writetocsv           - (logical) If true, writes intensity data to CSV.
%                              [Defaults to opts.writetocsv or false]
%       saveregisteredvolume - (logical) If true, saves the registered volume
%                              stack. [Defaults to opts.saveregisteredvol or false]
%       areafun              - (function_handle) Per-area summary statistic.
%                              [Defaults to opts.areafun or @median]
%
%   Outputs:
%       allmedians           - (single) Ngroups x 2 x Nchans array containing
%                              the per-area, per-hemisphere summary (the median
%                              unless 'areafun' says otherwise).
%
%   The per-channel .mat files keep the variable name 'medianoverareas' whatever
%   'areafun' is, so existing readers (loadMouseBackgroundSignal and friends)
%   keep working; the function actually used is saved alongside it as 'areafun'.
%--------------------------------------------------------------------------
% 1. Load configuration and transforms
%--------------------------------------------------------------------------
trstruct   = load(fullfile(savepath, 'transform_params.mat'));
loadedOpts = load(fullfile(savepath, 'regopts.mat'));
opts       = loadedOpts.opts;
%--------------------------------------------------------------------------
% 2. Parse Inputs and establish defaults
%--------------------------------------------------------------------------
% Determine defaults from the loaded opts structure, falling back to false
if isfield(opts, 'writetocsv')
    defaultWriteCsv = opts.writetocsv;
else
    defaultWriteCsv = false;
end

if isfield(opts, 'saveregisteredvol')
    defaultSaveVol  = opts.saveregisteredvol;
else
    defaultSaveVol  = false;
end

if isfield(opts, 'areafun') && isa(opts.areafun, 'function_handle')
    defaultAreaFun  = opts.areafun;
else
    defaultAreaFun  = @median;
end

% Set up the input parser
p = inputParser;
addRequired(p, 'savepath', @(x) ischar(x) || isstring(x));
addParameter(p, 'writetocsv', defaultWriteCsv, @(x) islogical(x) || isscalar(x));
addParameter(p, 'saveregisteredvolume', defaultSaveVol, @(x) islogical(x) || isscalar(x));
addParameter(p, 'areafun', defaultAreaFun, @(x) isa(x, 'function_handle'));

% Parse the arguments passed to the function
parse(p, savepath, varargin{:});

% Assign the final parsed variables for the rest of the script
writetocsv = p.Results.writetocsv;
saveregvol = p.Results.saveregisteredvolume;
areafun    = p.Results.areafun;

% accumarray needs a function that reduces a column vector to a scalar; check
% that here rather than letting it fail thousands of areas into the loop
try
    testout = areafun(single([1; 2; 3; 4]));
catch ME
    error('generateRegisteredBrainVolumes:badAreaFun', ...
        ['''areafun'' (%s) could not be applied to a column vector: %s. It must ' ...
         'take a vector of voxel values and return one number, like @median, ' ...
         '@mean or @std.'], func2str(areafun), ME.message);
end
if ~isscalar(testout) || ~isnumeric(testout)
    error('generateRegisteredBrainVolumes:badAreaFun', ...
        ['''areafun'' (%s) must return a single number per area; it returned a ' ...
         '%s of size %s.'], func2str(areafun), class(testout), mat2str(size(testout)));
end
%--------------------------------------------------------------------------
% 3. Initialize paths and channel names
%--------------------------------------------------------------------------
registerpath = fullfile(savepath, 'volume_registered');
makeNewDir(registerpath);

channames = cell(opts.Nchans, 1);
if isfield(opts, 'channames')
    channames = opts.channames;
end
%==========================================================================
fprintf('Applying transforms... \n'); savetic = tic;

straightvol = zeros([trstruct.atlassize opts.Nchans], 'uint16');

for ichan = 1:opts.Nchans
    %--------------------------------------------------------------------------
    % load and transform background volume
    volpath    = dir(fullfile(savepath, sprintf('chan_%d_*register*.tif', ichan)));
    isannot    = contains({volpath(:).name}, 'annotation');
    volpath(isannot) = [];
    currfname  = fullfile(volpath.folder, volpath.name);
    fprintf('Registering %s\n', currfname)
    %--------------------------------------------------------------------------
    backvolume = readDownStack(currfname);
    volume     = permuteBrainVolume(backvolume, trstruct.how_to_perm);

    volumereg  = transformix(volume,trstruct.tform_bspline_samp20um_to_atlas_20um_px,...
        'movingscale', opts.registres*1e-3*[1 1 1]);
    volumereg                   = uint16(abs(volumereg));
    Rmoving                     = imref3d(size(volumereg));
    Rfixed                      = imref3d(trstruct.atlassize);
    straightvol(:, :, :, ichan) = imwarp(volumereg, Rmoving, trstruct.tform_affine_samp20um_to_atlas_10um_px,...
        'OutputView',Rfixed);
    %--------------------------------------------------------------------------
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, opts.Nchans, toc(savetic));
end
%==========================================================================
if saveregvol
    fprintf('Saving registered volumes... '); savetic = tic;
    saveLargeSliceVolume(permute(straightvol, [1 2 4 3]), channames, registerpath);
    fprintf('Done! Took %2.2f s. \n', toc(savetic));
end
%==========================================================================
fprintf('Calculating background fluoresence in atlas coords (per-area %s)...\n', ...
    func2str(areafun)); proctic = tic;

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

Ngroups          = numel(areaidx);
Nforaccum        = max(av, [], 'all') + 1;
Npxlr            = size(av,3)/2;

allmedians = nan(Ngroups, 2, opts.Nchans, 'single');
for ichan = 1:opts.Nchans

    medianoverareas  = nan(Ngroups, 2, 'single');
    volumeoverareas  = nan(Ngroups, 2, 'single');
    for iside = 1:2
        istart = (iside - 1) * Npxlr + 1;
        iend   = istart + Npxlr - 1;
    
        sideav    = reshape(av(:, :, istart:iend), [], 1);
        sidevals  = reshape(straightvol(:, :, istart:iend, ichan), [], 1);
        ikeep     = sideav>0;
        medareas  = single(accumarray(sideav(ikeep)+1, sidevals(ikeep), [Nforaccum 1], areafun));
        medareas  = medareas(areaidx+1);

        volareas  = accumarray(sideav+1, 1, [Nforaccum 1], @sum);
        volareas  = volareas(areaidx+1);
        volumeoverareas(:, iside) = volareas * (opts.atlasres*1e-3)^3;

        % get index 0 level for background
        backlevel                 = single(areafun(sidevals(~ikeep)));
        medareas(areaidx == 0)    = backlevel;
        medianoverareas(:, iside) = medareas;
    end

    allmedians(:, :, ichan) = medianoverareas;

    % save as mat file for later processing
    fmatname  = fullfile(registerpath, sprintf('chan%02d_intensities.mat', ichan));
    areafunname = func2str(areafun);
    save(fmatname, 'medianoverareas', 'areaidx', 'volumeoverareas', 'areafun', 'areafunname')

    % save as csv if asked
    if writetocsv
        currtable = array2table([areaidx medianoverareas volumeoverareas], ...
            'VariableNames',...
            {'parcellation_index', ...
            'RightSideIntensity', 'LeftSideIntensity',...
             'RightSideVolume[mm3]', 'LeftSideVolume[mm3]'});
        currtable = addvars(currtable, namessub(ib), namesstruct(ibstr),  namesdiv(ibdiv),...
            'NewVariableNames',{'name', 'structure', 'division'},'Before','parcellation_index');
        fsavename      = fullfile(registerpath, sprintf('chan%02d_intensities.csv', ichan));
        writetable(currtable, fsavename, 'Delimiter', ';')
    end
    %--------------------------------------------------------------------------
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, opts.Nchans, toc(proctic));
    %--------------------------------------------------------------------------
end
%==========================================================================
end