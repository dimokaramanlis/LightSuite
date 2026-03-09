function allmedians = generateRegisteredBrainVolumes(savepath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
% we find useful files
trstruct   = load(fullfile(savepath, 'transform_params.mat'));
opts       = load(fullfile(savepath, 'regopts.mat'));
opts       = opts.opts;
writetocsv = getOr(opts, 'writetocsv', false);
saveregvol = getOr(opts, 'saveregisteredvol', false);
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
    currfname  = fullfile(volpath.folder, volpath.name);
    fprintf('Registering %s\n', currfname)
    %--------------------------------------------------------------------------
    backvolume = readDownStack(currfname);
    volume     = permute(backvolume, trstruct.how_to_perm);
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
fprintf('Calculating background fluoresence in atlas coords...\n'); proctic = tic;

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
parcelinfo       = readtable(fullfile(allen_atlas_path, 'parcellation_to_parcellation_term_membership.csv'));
substridx        = contains(parcelinfo.parcellation_term_set_name, 'substructure');
[areaidx, ib]    = unique(parcelinfo.parcellation_index(substridx));
namessub         = parcelinfo.parcellation_term_name(substridx);
Ngroups          = numel(areaidx);
Nforaccum        = max(av, [], 'all') + 1;
Npxlr            = size(av,3)/2;

allmedians = nan(Ngroups, 2, opts.Nchans, 'single');
for ichan = 1:opts.Nchans

    medianoverareas  = nan(Ngroups, 2, 'single');
    for iside = 1:2
        istart = (iside - 1) * Npxlr + 1;
        iend   = istart + Npxlr - 1;
    
        sideav    = reshape(av(:, :, istart:iend), [], 1);
        sidevals  = reshape(straightvol(:, :, istart:iend, ichan), [], 1);
        ikeep     = sideav>0;
        medareas  = single(accumarray(sideav(ikeep)+1, sidevals(ikeep), [Nforaccum 1], @median));
        medareas  = medareas(areaidx+1);
        % get index 0 level for background
        backlevel                 = single(median(sidevals(~ikeep)));
        medareas(areaidx == 0)    = backlevel;
        medianoverareas(:, iside) = medareas;
    end

    allmedians(:, :, ichan) = medianoverareas;

    % save as mat file for later processing
    fmatname  = fullfile(registerpath, sprintf('chan_%d_intensities.mat', ichan));
    save(fmatname, 'medianoverareas', 'areaidx')

    % save as csv if asked
    if writetocsv
        currtable = array2table([areaidx medianoverareas], ...
            'VariableNames',{'parcellation_index', 'RightSideIntensity', 'LeftSideIntensity'});
        currtable = addvars(currtable, namessub(ib), 'NewVariableNames','name','Before','parcellation_index');
        fsavename      = fullfile(registerpath, sprintf('chan_%d_intensities.csv', ichan));
        writetable(currtable, fsavename)
    end
    %--------------------------------------------------------------------------
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, opts.Nchans, toc(proctic));
    %--------------------------------------------------------------------------
end
%==========================================================================
end