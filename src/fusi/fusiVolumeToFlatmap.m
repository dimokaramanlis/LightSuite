function out = fusiVolumeToFlatmap(vol, opts)
%FUSIVOLUMETOFLATMAP Project a fUSI atlas-space volume onto an Allen flatmap.
%   out = fusiVolumeToFlatmap(vol, opts) takes a single volume in Allen CCF
%   space (e.g. one of the atlas-space activation / VOSI maps produced by
%   fusiMouseActivationMaps) and returns its whole-cortex "flatbrain"
%   projection - the average of the volume along each cortical streamline -
%   using the Allen ccf_streamlines / allensdk Python packages.
%
%   MATLAB does not do the projection itself: it hands the volume to the
%   Python worker src/fusi/py/fusi_flatmap_project.py, run with the
%   interpreter given in opts.pythonexe (the conda/venv env that has allensdk
%   and ccf_streamlines installed), and reads the 2-D result back. Data is
%   exchanged as small v7 .mat files, so the interpreter can live in any
%   environment (or machine) without linking into MATLAB.
%
%   The volume must be in CCF orientation (dim1=AP, dim2=DV, dim3=ML),
%   isotropic at opts.inputres_um. A 50 um volume is therefore [264 160 228],
%   which is exactly size(optA.atlas) for the fUSI atlas (see
%   prepare_fusi_atlas.m); the maps in <mouse>_mouse_actmaps.mat (atlasMaps.*,
%   vosiMaps.*) are already in this space and can be projected directly. The
%   worker upsamples the volume to the 10 um grid the streamlines are defined
%   on, so no resampling is needed on the MATLAB side.
%
%   Input:
%     vol   3-D numeric volume in CCF (AP x DV x ML) at opts.inputres_um.
%     opts:
%        .pythonexe   REQUIRED full path to the env's python(.exe) that has
%                     allensdk + ccf_streamlines.
%        .resourcedir dir holding the shared flatmap resources
%                     (flatmap_butterfly.h5 / .nrrd, surface_paths_10_v3.h5,
%                     labelDescription_ITKSNAPColor.txt, manifest.json).
%                     Default 'D:\AllenAtlas\flatmaps'.
%        .inputres_um isotropic voxel size of vol in um (default 50).
%        .hemisphere  'left' (default) or 'right'.
%        .kind        streamline reduction: 'mean' (default) or 'max'.
%        .maskres_um  annotation resolution for the cortex mask (default 50).
%        .boundaries  also return the area-boundary overlay image (default true).
%        .scriptpath  path to the Python worker (default: py/ next to this file).
%        .tempdir     scratch dir for the exchange files (default tempdir).
%        .keeptemp    keep the exchange .mat files (default false).
%        .verbose     default true.
%
%   Output struct out:
%     .flatmap     2-D single flatmap (streamline-averaged cortex).
%     .boundaries  2-D uint8 area-boundary image (same size), or [].
%     .regionBoundaries struct with .names {1xN acronym} and .coords {1xN Nx2
%                  [row col] flatmap pixels}, so a chosen subset of areas can be
%                  outlined (e.g. the dorsal/ventral visual streams).
%     .hemisphere, .kind, .inputres_um, .maskres_um   settings used.
%
%   Example:
%     r = load('DS_WT61_mouse_actmaps.mat');
%     o = struct('pythonexe', 'C:\Users\me\miniconda3\envs\ccf\python.exe');
%     fm = fusiVolumeToFlatmap(r.vosiMaps.full, o);
%     imagesc(fm.flatmap.'); axis image off; colormap magma
%
%   See also FUSIMOUSEACTIVATIONMAPS, PREPARE_FUSI_ATLAS.

if nargin < 2, opts = struct(); end
opts = setdefault(opts, 'pythonexe',  '');
opts = setdefault(opts, 'resourcedir', 'D:\AllenAtlas\flatmaps');
opts = setdefault(opts, 'inputres_um', 50);
opts = setdefault(opts, 'hemisphere', 'left');
opts = setdefault(opts, 'kind',       'mean');
opts = setdefault(opts, 'maskres_um', 50);
opts = setdefault(opts, 'boundaries', true);
opts = setdefault(opts, 'scriptpath', fullfile(fileparts(mfilename('fullpath')), ...
    'py', 'fusi_flatmap_project.py'));
opts = setdefault(opts, 'tempdir',    tempdir);
opts = setdefault(opts, 'keeptemp',   false);
opts = setdefault(opts, 'verbose',    true);

% ---- validate inputs ---------------------------------------------------
if isempty(opts.pythonexe)
    error('fusiVolumeToFlatmap:noPython', ...
        ['opts.pythonexe is required: give the full path to the python(.exe) ', ...
         'of the environment that has allensdk + ccf_streamlines.']);
end
if exist(opts.pythonexe, 'file') ~= 2
    error('fusiVolumeToFlatmap:badPython', 'python not found: %s', opts.pythonexe);
end
if exist(opts.scriptpath, 'file') ~= 2
    error('fusiVolumeToFlatmap:noScript', 'worker script not found: %s', opts.scriptpath);
end
if exist(opts.resourcedir, 'dir') ~= 7
    error('fusiVolumeToFlatmap:noResources', 'resourcedir not found: %s', opts.resourcedir);
end
if ndims(vol) ~= 3
    error('fusiVolumeToFlatmap:notVolume', 'vol must be a 3-D volume.');
end

% shape sanity check: a res-um CCF (AP,DV,ML) volume has this many voxels.
expshape = round([1320 800 1140] * 10 / opts.inputres_um);
if ~isequal(size(vol, 1:3), expshape)
    error('fusiVolumeToFlatmap:shape', ...
        ['vol is %s but a %g um CCF (AP,DV,ML) volume should be %s. Check the ', ...
         'resolution/orientation (dim1=AP, dim2=DV, dim3=ML).'], ...
        mat2str(size(vol, 1:3)), opts.inputres_um, mat2str(expshape));
end

vol         = single(vol);            %#ok<NASGU> saved below
inputres_um = double(opts.inputres_um); %#ok<NASGU>

% ---- exchange files ----------------------------------------------------
tag     = sprintf('%s_%09d', datestr(now, 'yyyymmddHHMMSS'), randi(1e9)); %#ok<TNOW1,DATST>
infile  = fullfile(opts.tempdir, sprintf('fusiflat_in_%s.mat',  tag));
outfile = fullfile(opts.tempdir, sprintf('fusiflat_out_%s.mat', tag));
cleaner = onCleanup(@() cleanupTemp(infile, outfile, opts.keeptemp));

save(infile, 'vol', 'inputres_um', '-v7');   % v7 so scipy.io can read it

% ---- build + run the command ------------------------------------------
args = sprintf(['"%s" --input "%s" --output "%s" --resourcedir "%s" ', ...
    '--hemisphere %s --kind %s --maskres %d'], ...
    opts.scriptpath, infile, outfile, opts.resourcedir, ...
    opts.hemisphere, opts.kind, round(opts.maskres_um));
if ~opts.boundaries, args = [args ' --no-boundaries']; end

cmd = sprintf('"%s" %s', opts.pythonexe, args);
if ispc
    % Prepend the env dirs to PATH so conda DLLs (e.g. VTK used by
    % ccf_streamlines) resolve when calling python.exe directly.
    envroot = fileparts(opts.pythonexe);
    pathadd = strjoin({envroot, fullfile(envroot, 'Library', 'bin'), ...
        fullfile(envroot, 'Library', 'mingw-w64', 'bin'), ...
        fullfile(envroot, 'Library', 'usr', 'bin'), ...
        fullfile(envroot, 'Scripts'), fullfile(envroot, 'bin')}, ';');
    cmd = sprintf('set "PATH=%s;%%PATH%%" && %s', pathadd, cmd);
end

if opts.verbose, fprintf('fusiVolumeToFlatmap: running Python worker...\n'); end
[status, cmdout] = system(cmd);
if opts.verbose && ~isempty(strtrim(cmdout)), fprintf('%s\n', strtrim(cmdout)); end
if status ~= 0
    error('fusiVolumeToFlatmap:pythonFailed', ...
        'Python worker failed (status %d):\n%s', status, cmdout);
end
if exist(outfile, 'file') ~= 2
    error('fusiVolumeToFlatmap:noOutput', ...
        'Python worker produced no output. Log:\n%s', cmdout);
end

% ---- collect results ---------------------------------------------------
R = load(outfile);
out = struct();
out.flatmap     = single(R.flatmap);
out.boundaries  = [];
if isfield(R, 'boundaries'), out.boundaries = R.boundaries; end
% per-region boundary outlines (acronym -> Nx2 [row col] flatmap pixels), so a
% chosen subset of areas can be highlighted rather than only the merged overlay.
out.regionBoundaries = struct('names', {{}}, 'coords', {{}});
if isfield(R, 'region_names') && isfield(R, 'region_coords')
    nm = R.region_names; co = R.region_coords;
    if ~iscell(nm), nm = num2cell(nm); end
    if ~iscell(co), co = num2cell(co); end
    out.regionBoundaries.names  = cellfun(@(s) char(string(s)), nm(:).', 'UniformOutput', false);
    out.regionBoundaries.coords = cellfun(@double, co(:).', 'UniformOutput', false);
end
out.hemisphere  = opts.hemisphere;
out.kind        = opts.kind;
out.inputres_um = opts.inputres_um;
out.maskres_um  = opts.maskres_um;
if opts.verbose
    fprintf('fusiVolumeToFlatmap: flatmap %s (%s, %s).\n', ...
        mat2str(size(out.flatmap)), opts.hemisphere, opts.kind);
end
end

% =========================================================================
function cleanupTemp(infile, outfile, keep)
if keep, return; end
if exist(infile,  'file') == 2, delete(infile);  end
if exist(outfile, 'file') == 2, delete(outfile); end
end

function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
