function [cf, cf_intensity, pp, resout] = plotGRINTracingResults(savepath, varargin)
% PLOTGRINTRACINGRESULTS  Load and display GRIN lens tracing results saved
%   by annotateGRINLens in a LightSuite folder.
%
%   [CF, PP, RESOUT] = plotGRINTracingResults(SAVEPATH)
%
%   Searches SAVEPATH for every  grin_fiber<N>_atlas.mat  file produced by
%   annotateGRINLens (one per annotated fiber), loads them in ascending
%   fiber-number order, then:
%
%   1. Re-samples the atlas cross-sections at the standard depth range
%      (grinTracingDepths, -400:100:400 µm) for any fiber saved at a narrower
%      range, so the negative offsets above the tip — used to check that the
%      fit follows the visible fiber track — are present too.
%
%   2. Calls plotGRINAtlasImages to show:
%       • a 3-D brain-grid plot with colour-coded fiber cylinders, and
%       • one row of circular parcellation-region panels per fiber at the
%         non-negative depths (0–400 µm). The negative depths are stored but
%         not drawn.
%
%   3. Produces a companion figure (plotGRINIntensityImages) of the
%      fluorescence cross-sections for every channel at the same depths.
%      By default ('Resolution','full') they are sampled from the RAW
%      full-resolution data (extractGRINRawIntensitySlices), so single cells
%      are visible. With 'Resolution','registered', or when the raw data
%      cannot be found, they come from SAVEPATH/volume_registered/ (20 µm data
%      on the 10 µm atlas grid, extractGRINIntensitySlices).
%
%      Full-resolution slices are stored in the fiber files and reused on
%      later calls as long as the fit, depths, pixel size and channels are
%      unchanged, so the raw data is read only once.
%
%   4. Writes the regenerated depths and the intensity cross-sections back
%      into each grin_fiber<N>_atlas.mat, so both live alongside the fit:
%
%        depths_um         1×Ndepths depth offsets (µm), now including negatives
%        slices_av         1×Ndepths cell of annotation cross-sections
%        rvec_arr          1×Ndepths cell of coordinate vectors
%        slices_int        1×Ndepths cell, each [Ngrid×Ngrid×Nchan] single
%        rvec_int          1×Ndepths cell, grid of slices_int (atlas voxels)
%        slices_av_int     1×Ndepths cell, annotation on the slices_int grid
%                          (full resolution only)
%        median_intensity  [Ndepths × Nchan], median inside the implant disc
%        mean_intensity    [Ndepths × Nchan], mean inside the implant disc
%        chan_names        1×Nchan channel names, labelling the 3rd dimension
%        int_meta          what the slices were sampled from (cache key)
%
%      The file is rewritten in place (atomically) rather than appended to,
%      the fit is left untouched, and any other variables are preserved. A
%      file is only rewritten when its depths were regenerated or intensity
%      was added.
%
%   The Allen CCF atlas resources (annotation_10.nii.gz and
%   parcellation_to_parcellation_term_membership.csv) must be on the
%   MATLAB path, exactly as required by annotateGRINLens.
%
%   Usage
%   -----
%       plotGRINTracingResults(savepath)
%       plotGRINTracingResults(savepath, 'Resolution', 'registered')
%       plotGRINTracingResults(savepath, 'PixelSize', 2, 'Channels', 2)
%       [cf, cf_intensity, pp, resout] = plotGRINTracingResults(savepath)
%
%   Inputs
%   ------
%   savepath  – LightSuite output folder in which annotateGRINLens was run.
%               Must contain at least one  grin_fiber*_atlas.mat  file.
%
%   Name-value options
%   ------------------
%   'Resolution' ('full')  'full' = raw data, 'registered' = volume_registered/.
%   'PixelSize'  ([])      full-res grid spacing in µm (default: raw pixel size).
%   'Channels'   ([])      raw channels to sample (default: all).
%   'DataFolder' ('')      raw data folder, if it moved since processing.
%   'Recompute'  (false)   resample even if stored slices match.
%
%   Outputs (optional)
%   ------------------
%   cf           – handle to the atlas-annotation figure (plotGRINAtlasImages).
%   cf_intensity – handle to the intensity figure ([] if volume_registered/
%                  is absent or intensity extraction failed).
%   pp           – panel object returned by plotGRINAtlasImages.
%   resout – struct with fields:
%               .all_results        – loaded fiber atlas structs
%               .intensity_results  – per-fiber intensity structs (empty if
%                                     volume_registered/ is absent)
%               .channames          – channel name strings (empty if above)
%               .cf_atlas           – atlas figure handle
%               .cf_intensity       – intensity figure handle ([] if absent)
%
%   Example
%   -------
%       plotGRINTracingResults('/data/mouse01/lightsuite_output');

    ip = inputParser;
    addParameter(ip, 'Resolution', 'full', @(x) any(strcmpi(x, {'full', 'registered'})));
    addParameter(ip, 'PixelSize',  [],     @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(ip, 'Channels',   [],     @isnumeric);
    addParameter(ip, 'DataFolder', '',     @(x) ischar(x) || isstring(x));
    addParameter(ip, 'Recompute',  false,  @(x) islogical(x) || isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    %----------------------------------------------------------------------
    % 1. Find saved per-fiber result files
    %----------------------------------------------------------------------
    pattern   = fullfile(savepath, 'grin_fiber*_atlas.mat');
    flist     = dir(pattern);

    if isempty(flist)
        error('plotGRINTracingResults: no grin_fiber*_atlas.mat files found in:\n  %s', savepath);
    end

    % Sort by the embedded fiber number (e.g. grin_fiber3_atlas.mat → 3)
    fiber_nums = zeros(1, numel(flist));
    for k = 1:numel(flist)
        tok = regexp(flist(k).name, 'grin_fiber(\d+)_atlas', 'tokens', 'once');
        if ~isempty(tok)
            fiber_nums(k) = str2double(tok{1});
        end
    end
    [~, sort_idx] = sort(fiber_nums);
    flist = flist(sort_idx);

    fprintf('plotGRINTracingResults: found %d fiber result file(s) in\n  %s\n', ...
        numel(flist), savepath);

    %----------------------------------------------------------------------
    % 2. Load per-fiber structs
    %----------------------------------------------------------------------
    all_results = cell(1, numel(flist));
    fpaths      = cell(1, numel(flist));
    for k = 1:numel(flist)
        fpath      = fullfile(flist(k).folder, flist(k).name);
        fpaths{k}  = fpath;
        fprintf('  Loading %s ...\n', flist(k).name);
        s = load(fpath);

        % Validate that the required fields are present
        required = {'center_vox', 'normal_vox', 'radius_vox', ...
                    'atlas_pts', 'depths_um', 'slices_av', 'rvec_arr'};
        missing = required(~isfield(s, required));
        if ~isempty(missing)
            error('plotGRINTracingResults: file %s is missing field(s): %s', ...
                flist(k).name, strjoin(missing, ', '));
        end

        all_results{k} = s;
    end

    %----------------------------------------------------------------------
    % 3. Load Allen CCF atlas resources (annotation volume + parcellation)
    %----------------------------------------------------------------------
    allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
    if isempty(allen_atlas_path)
        error('plotGRINTracingResults: annotation_10.nii.gz not found on the MATLAB path.');
    end

    parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));
    if isempty(parcel_path)
        error('plotGRINTracingResults: parcellation_to_parcellation_term_membership.csv not found on the MATLAB path.');
    end

    fprintf('Loading parcellation table ...\n');
    parcelinfo = readtable(fullfile(parcel_path, 'parcellation_to_parcellation_term_membership.csv'));

    substridx   = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib] = unique(parcelinfo.parcellation_index(substridx));
    namessub    = parcelinfo.parcellation_term_acronym(substridx);
    namessub    = namessub(ib);
    colorsub    = [parcelinfo.red(substridx), ...
                   parcelinfo.green(substridx), ...
                   parcelinfo.blue(substridx)];
    colorsub    = colorsub(ib, :);

    %----------------------------------------------------------------------
    % 4. Bring every fiber onto the standard depth range
    %----------------------------------------------------------------------
    % Files saved by an older annotateGRINLens hold only the non-negative
    % depths. Re-sample the atlas cross-sections at the full range (including
    % the negative offsets above the tip) so the stored data is uniform; this
    % needs the annotation volume, but only when something is actually missing.
    % Files saved with the older, narrower slice extent are re-sampled too.
    av_file = fullfile(allen_atlas_path, 'annotation_10.nii.gz');
    [all_results, atlas_regen, av] = ensureStandardDepths(all_results, av_file);

    %----------------------------------------------------------------------
    % 5. Intensity: raw full-resolution data, or volume_registered/
    %----------------------------------------------------------------------
    cf_intensity       = [];
    intensity_results  = {};
    channames_int      = {};
    intensity_added    = false(1, numel(all_results));
    int_source         = '';

    if strcmpi(opt.Resolution, 'full')
        % pixel size / channels the raw sampler would use, for the cache key
        ro   = load(fullfile(savepath, 'regopts.mat'), 'opts');
        tp   = load(fullfile(savepath, 'transform_params.mat'), 'ori_pxsize');
        pxum = opt.PixelSize;
        if isempty(pxum); pxum = min(tp.ori_pxsize); end
        chans = opt.Channels;
        if isempty(chans); chans = 1:getOr(ro.opts, 'Nchans', 1); end

        metas = cell(1, numel(all_results));
        for k = 1:numel(all_results)
            r = all_results{k};
            metas{k} = struct('source', 'raw', 'pxsize_um', pxum, ...
                'channels', chans(:)', 'center_vox', r.center_vox(:)', ...
                'normal_vox', r.normal_vox(:)', 'radius_vox', r.radius_vox, ...
                'depths_um', r.depths_um(:)', 'grid', grinSliceGrid(r.radius_vox, pxum / 10));
        end
        cached = ~opt.Recompute && all(cellfun(@(r, m) isfield(r, 'int_meta') && ...
            isfield(r, 'rvec_int') && isfield(r, 'slices_av_int') && ...
            isequal(r.int_meta, m), all_results, metas));

        if cached
            fprintf('Using stored full-resolution intensity slices (%.2f µm).\n', pxum);
            intensity_results = cell(1, numel(all_results));
            for k = 1:numel(all_results)
                r = all_results{k};
                intensity_results{k} = struct('slices_int', {r.slices_int}, ...
                    'median_intensity', r.median_intensity, ...
                    'mean_intensity', r.mean_intensity, 'rvec_arr', {r.rvec_int}, ...
                    'depths_um', r.depths_um, 'radius_vox', r.radius_vox, ...
                    'pxsize_um', pxum, 'slices_av', {r.slices_av_int});
            end
            channames_int = all_results{1}.chan_names;
            int_source    = 'raw';
        else
            fprintf('Extracting full-resolution intensity slices from the raw data ...\n');
            try
                if isempty(av)   % annotation on the same fine grid needs the volume
                    fprintf('Loading annotation volume ...\n');
                    av = niftiread(av_file);
                end
                [intensity_results, channames_int] = extractGRINRawIntensitySlices( ...
                    savepath, all_results, 'PixelSize', pxum, 'Channels', chans, ...
                    'DataFolder', opt.DataFolder, 'Annotation', av);
                int_source      = 'raw';
                intensity_added = true(1, numel(all_results));
                for k = 1:numel(all_results)
                    all_results{k}.int_meta = metas{k};
                end
            catch ME
                warning('plotGRINTracingResults:rawFailed', ...
                    'full-resolution extraction failed, using volume_registered/:\n  %s', ...
                    ME.message);
            end
        end
    end

    if isempty(int_source)
        if exist(fullfile(savepath, 'volume_registered'), 'dir')
            fprintf('Found volume_registered/ – extracting intensity profiles ...\n');
            try
                [intensity_results, channames_int] = ...
                    extractGRINIntensitySlices(savepath, all_results);
                int_source      = 'registered';
                intensity_added = true(1, numel(all_results));
                for k = 1:numel(all_results)
                    all_results{k}.int_meta = struct('source', 'registered');
                end
            catch ME
                warning('plotGRINTracingResults:intensityFailed', ...
                    'intensity extraction failed:\n  %s', ME.message);
            end
        else
            fprintf('(volume_registered/ not found – skipping intensity figure)\n');
        end
    end

    clear av;

    % Carry the sampled images on the results so they are stored with the
    % fit (step 7) and returned to the caller.
    for k = find(intensity_added)
        all_results{k}.slices_int       = intensity_results{k}.slices_int;
        all_results{k}.rvec_int         = intensity_results{k}.rvec_arr;
        all_results{k}.median_intensity = intensity_results{k}.median_intensity;
        all_results{k}.mean_intensity   = intensity_results{k}.mean_intensity;
        all_results{k}.chan_names       = channames_int;
        if strcmp(int_source, 'raw')
            all_results{k}.slices_av_int = intensity_results{k}.slices_av;
        end
    end

    %----------------------------------------------------------------------
    % 6. Figures (non-negative depths only)
    %----------------------------------------------------------------------
    % At full resolution the atlas figure shows the annotation on the same
    % fine grid as the intensity; otherwise the stored 10 µm slices.
    plot_results = all_results;
    if strcmp(int_source, 'raw')
        for k = 1:numel(plot_results)
            plot_results{k}.slices_av = intensity_results{k}.slices_av;
            plot_results{k}.rvec_arr  = intensity_results{k}.rvec_arr;
        end
    end

    fprintf('Plotting atlas annotation figure (%d fiber(s)) ...\n', numel(all_results));
    [cf, pp, resout] = plotGRINAtlasImages(plot_results, areaidx, namessub, colorsub);
    if ~isempty(cf)
        cf.Name = sprintf('GRIN Atlas Regions – %s', savepath);
    end

    if ~isempty(int_source)
        [cf_intensity, ~] = plotGRINIntensityImages( ...
            plot_results, intensity_results, channames_int);
        if ~isempty(cf_intensity)
            cf_intensity.Name = sprintf('GRIN Intensity Profiles (%s) – %s', ...
                int_source, savepath);
        end
    end

    %----------------------------------------------------------------------
    % 7. Persist the regenerated depths and the intensity images
    %----------------------------------------------------------------------
    % Written whenever the atlas depths were regenerated or intensity was
    % added, so the negative-depth cross-sections and the sampled fluorescence
    % live alongside the fit they came from.
    storeFiberResults(fpaths, all_results, atlas_regen | intensity_added);

    %----------------------------------------------------------------------
    % 8. Pack optional output struct
    %----------------------------------------------------------------------
    if nargout >= 3
        resout.all_results       = all_results;
        resout.intensity_results = intensity_results;
        resout.channames         = channames_int;
        resout.cf_atlas          = cf;
        resout.cf_intensity      = cf_intensity;
    end
end

%==========================================================================

function [all_results, regen, av] = ensureStandardDepths(all_results, av_file)
% Re-sample each fiber's atlas cross-sections at grinTracingDepths() and on
% the grinSliceGrid() extent when the stored ones differ, and report which
% fibers were changed.
%
% The annotation volume is loaded once, and only if at least one fiber needs
% regenerating, so folders already up to date pay nothing. It is returned
% (empty if not loaded) so the caller can reuse it.
    standard = grinTracingDepths();
    ns       = numel(standard);
    n        = numel(all_results);
    regen    = false(1, n);
    av       = [];

    for k = 1:n
        r = all_results{k};
        d = r.depths_um(:)';
        if ~isequal(d, standard) || numel(r.slices_av) ~= ns || ...
                ~isequal(r.rvec_arr{1}(:)', grinSliceGrid(r.radius_vox, 1))
            regen(k) = true;
        end
    end
    if ~any(regen); return; end

    fprintf('Re-sampling atlas cross-sections at %d depths (%d to %d µm) for %d fiber(s) ...\n', ...
        ns, standard(1), standard(end), nnz(regen));
    av = niftiread(av_file);
    dvox = standard / 10;

    for k = find(regen)
        r  = all_results{k};
        sl = cell(1, ns);
        rv = cell(1, ns);
        for ii = 1:ns
            [sl{ii}, rv{ii}] = extractAnnotationCircularSlice( ...
                av, r.center_vox, r.normal_vox, r.radius_vox, dvox(ii));
        end
        r.depths_um = standard;
        r.slices_av = sl;
        r.rvec_arr  = rv;
        all_results{k} = r;
    end
end

function storeFiberResults(fpaths, all_results, dirty)
% Rewrite each changed grin_fiber<N>_atlas.mat with the regenerated depth
% range and, when present, the sampled intensity images.
%
% The file is read, updated and rewritten rather than appended to. MATLAB's
% '-append' does not reclaim the records it replaces in a compressed v7
% MAT-file, so re-running this function on the same folder grew the file by
% roughly its own size every time. Everything already in the file is loaded
% first, so unknown fields survive the round trip and the fit itself is only
% touched where it is meant to be (depths_um, slices_av, rvec_arr).
%
% The rewrite goes to a temporary file that is then moved into place, so an
% interrupted save cannot leave a half-written result file where a complete
% one used to be. A file that cannot be written (read-only, open elsewhere,
% on a disconnected share) is reported and skipped: the figures are already
% drawn by this point, and losing them to a save error would be worse than
% not storing the extra data.
    fields = {'depths_um', 'slices_av', 'rvec_arr', ...
              'slices_int', 'rvec_int', 'slices_av_int', 'median_intensity', 'mean_intensity', ...
              'chan_names', 'int_meta'};

    for k = 1:numel(fpaths)
        if ~dirty(k); continue; end

        tmp = [fpaths{k} '.saving'];
        try
            S = load(fpaths{k});
            for f = fields
                if isfield(all_results{k}, f{1})
                    S.(f{1}) = all_results{k}.(f{1});
                end
            end
            save(tmp, '-struct', 'S');
            movefile(tmp, fpaths{k}, 'f');
        catch ME
            if exist(tmp, 'file'); delete(tmp); end
            warning('plotGRINTracingResults:resultsNotStored', ...
                'could not update %s:\n  %s', fpaths{k}, ME.message);
            continue
        end

        [~, nm, ex] = fileparts(fpaths{k});
        if isfield(all_results{k}, 'slices_int')
            extra = sprintf(', %d-channel intensity', numel(all_results{k}.chan_names));
        else
            extra = '';
        end
        fprintf('  Updated %s%s (%d depths stored%s)\n', ...
            nm, ex, numel(all_results{k}.depths_um), extra);
    end
end
