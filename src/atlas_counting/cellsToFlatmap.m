function out = cellsToFlatmap(cellcoords, opts)
%CELLSTOFLATMAP Map atlas-space cell detections onto the Allen cortical flatmap.
%   out = cellsToFlatmap(cellcoords, opts) gives every cell its own
%   continuous flatmap position, its depth below the pia along the cortical
%   streamline and its cortical layer, keeping the two hemispheres separate on
%   the two wings of the butterfly flatmap. Cells are NOT binned into a volume and
%   averaged along streamlines (as fusiVolumeToFlatmap does for volumes), so
%   no resolution is lost to the projection itself; the 10 um streamline grid
%   and the registration error set the actual precision.
%
%   How: for every 10 um voxel crossed by the Allen cortical streamlines
%   (ccf_streamlines surface_paths_10_v3.h5) a compact table stores the
%   flatmap position of a representative streamline, the voxel's depth along
%   it and the streamline id (for its layer extents). A cell takes the
%   trilinear interpolation of position and depth over its 8 surrounding
%   voxels, or - where streamlines leave gaps - a Gaussian-weighted average
%   over its 3x3x3 neighbourhood (highest-weight voxel where the flattening
%   tears). Right-hemisphere cells are looked up in the mirrored table and
%   placed on the right wing.
%
%   The work runs in Python (worker src/atlas_counting/py/cells_flatmap_project.py)
%   with the interpreter in opts.pythonexe, exchanging v7 .mat files, as in
%   fusiVolumeToFlatmap. The first call builds the voxel table once
%   (~1-2 min, ~0.7 GB on disk in <resourcedir>/cells_flatmap_cache); later
%   calls memory-map it and process cells in batches, e.g. ~6 s and ~0.6 GB
%   of private memory for 545k cells.
%
%   Input:
%     cellcoords  N x >=3 cell coordinates. Extra columns are ignored.
%                 opts.coordtype 'atlaspx' (default): LightSuite atlas pixels
%                 as saved in cell_locations_atlas.mat (atlasptcoords), i.e.
%                 [DV AP ML] = [col row plane] of the 10 um annotation volume,
%                 1-based, pixel centres at integers.
%                 opts.coordtype 'ccf_um': Allen CCF microns [AP DV ML].
%     opts:
%        .pythonexe   REQUIRED path to the python(.exe) with ccf_streamlines.
%        .resourcedir flatmap resources (flatmap_butterfly.h5/.nrrd,
%                     surface_paths_10_v3.h5, cortical_layers_10_v2.h5,
%                     avg_layer_depths.json, labelDescription_ITKSNAPColor.txt).
%                     Default 'D:\AllenAtlas\flatmaps'.
%        .coordtype   'atlaspx' (default) or 'ccf_um'.
%        .atlasres_um voxel size of 'atlaspx' coordinates (default 10).
%        .tearpx      flatmap spread (px) of the contributing voxels above which
%                     the highest-weight voxel is used instead of interpolating
%                     (default 10).
%        .minweight   trilinear weight the streamline voxels must carry to use
%                     trilinear interpolation (default 0.5); below it the 3x3x3
%                     Gaussian fallback is used.
%        .sigma       width of that fallback in voxels (default 0.75).
%        .batch       cells per batch in the worker (default 100000; lower
%                     it to reduce memory further).
%        .layermaps   also return flatmaps of local layer depths (default false).
%        .boundaries  also return area outlines (default true).
%        .scriptpath  worker path (default: py/ next to this file).
%        .tempdir     scratch dir for exchange files (default tempdir).
%        .keeptemp    keep exchange files (default false).
%        .verbose     default true.
%
%   Output struct out (N = number of cells, NaN / 0 for unmapped cells):
%     .flatxy        N x 2 butterfly-flatmap position [x y] in MATLAB pixel
%                    units, so an image F of size(.viewmask) shown with
%                    imagesc(F.') lines up with plot(flatxy(:,1), flatxy(:,2)).
%                    Left-hemisphere cells have x < .midline_x, right ones
%                    x > .midline_x (mirror image). The left wing is the frame
%                    of fusiVolumeToFlatmap (hemisphere 'left'). 1 px = 10 um.
%     .depth_um      N x 1 depth below pia along the streamline (um).
%     .depthnorm_um  N x 1 depth in the reference cortex, where every layer has
%                    its average thickness (comparable across areas).
%     .layer         N x 1 layer index into .layers.names (0 = unmapped).
%     .layeredges_um N x 7 local layer edges along the cell's own streamline
%                    (um below pia): tops of L1 L2/3 L4 L5 L6a L6b, then the
%                    white-matter border. NaN where a layer is absent (L4 in
%                    agranular areas such as PL).
%     .layers        reference layers: .names {1x6}, .top_um, .bottom_um
%                    (the frame of depthnorm_um).
%     .mapped        N x 1 logical, cell received a flatmap position. Cells
%                    with no streamline voxel within one voxel stay unmapped
%                    (NaN); cells just outside isocortex (white matter, pia)
%                    can still be mapped, so filter by annotation if needed.
%     .isright       N x 1 logical, cell is in the right hemisphere (CCF ML
%                    beyond the midline).
%     .tear          N x 1 logical, nearest voxel used (flattening tear).
%     .fallback      N x 1 logical, 3x3x3 fallback used instead of trilinear.
%     .viewmask      butterfly image, true where cortex exists (display with .').
%     .midline_x     x separating the two wings.
%     .layermaps_um  (if opts.layermaps) size(viewmask) x 7 local layer edges
%                    per flatmap pixel, as in .layeredges_um.
%     .regionBoundaries  .names {1xR}, .coords {1xR Mx2 [x y]}, .isright (1xR)
%                    area outlines of both wings, in the units of flatxy.
%     .coordtype     setting used.
%
%   Example:
%     s  = load('D:\DATA_folder\Mice\YX003\Anatomy\cell_locations_atlas.mat');
%     o  = struct('pythonexe', 'C:\...\envs\allensdk\python.exe');
%     fm = cellsToFlatmap(s.atlasptcoords, o);
%     plot(fm.flatxy(:,1), fm.flatxy(:,2), '.', 'MarkerSize', 1); axis ij image
%   See demos/ls_cells_flatmap.m for a full example.
%
%   See also FUSIVOLUMETOFLATMAP.

if nargin < 2, opts = struct(); end
opts = setdefault(opts, 'pythonexe',   '');
opts = setdefault(opts, 'resourcedir', 'D:\AllenAtlas\flatmaps');
opts = setdefault(opts, 'coordtype',   'atlaspx');
opts = setdefault(opts, 'atlasres_um', 10);
opts = setdefault(opts, 'tearpx',      10);
opts = setdefault(opts, 'minweight',   0.5);
opts = setdefault(opts, 'sigma',       0.75);
opts = setdefault(opts, 'batch',       100000);
opts = setdefault(opts, 'layermaps',   false);
opts = setdefault(opts, 'boundaries',  true);
opts = setdefault(opts, 'scriptpath',  fullfile(fileparts(mfilename('fullpath')), ...
    'py', 'cells_flatmap_project.py'));
opts = setdefault(opts, 'tempdir',     tempdir);
opts = setdefault(opts, 'keeptemp',    false);
opts = setdefault(opts, 'verbose',     true);

% ---- validate inputs ---------------------------------------------------
if isempty(opts.pythonexe)
    error('cellsToFlatmap:noPython', ...
        'opts.pythonexe is required: the python(.exe) of the env with ccf_streamlines.');
end
if exist(opts.pythonexe, 'file') ~= 2
    error('cellsToFlatmap:badPython', 'python not found: %s', opts.pythonexe);
end
if exist(opts.scriptpath, 'file') ~= 2
    error('cellsToFlatmap:noScript', 'worker script not found: %s', opts.scriptpath);
end
if exist(opts.resourcedir, 'dir') ~= 7
    error('cellsToFlatmap:noResources', 'resourcedir not found: %s', opts.resourcedir);
end
if ~ismatrix(cellcoords) || size(cellcoords, 2) < 3
    error('cellsToFlatmap:badCoords', 'cellcoords must be N x 3 (or wider).');
end

% ---- to CCF microns [AP DV ML] -----------------------------------------
xyz = double(cellcoords(:, 1:3));
switch lower(opts.coordtype)
    case 'atlaspx'
        % pixel p (1-based, centre at p) spans [(p-1) p) voxels -> centre (p-0.5)*res
        coords_um = single((xyz(:, [2 1 3]) - 0.5) * opts.atlasres_um); %#ok<NASGU>
    case 'ccf_um'
        coords_um = single(xyz); %#ok<NASGU>
    otherwise
        error('cellsToFlatmap:coordtype', 'unknown coordtype ''%s''.', opts.coordtype);
end

% ---- exchange files ----------------------------------------------------
tag     = sprintf('%s_%09d', datestr(now, 'yyyymmddHHMMSS'), randi(1e9)); %#ok<TNOW1,DATST>
infile  = fullfile(opts.tempdir, sprintf('cellsflat_in_%s.mat',  tag));
outfile = fullfile(opts.tempdir, sprintf('cellsflat_out_%s.mat', tag));
cleaner = onCleanup(@() cleanupTemp(infile, outfile, opts.keeptemp));
save(infile, 'coords_um', '-v7');

args = sprintf(['"%s" --input "%s" --output "%s" --resourcedir "%s" ', ...
    '--tearpx %g --minweight %g --sigma %g --batch %d'], ...
    opts.scriptpath, infile, outfile, opts.resourcedir, ...
    opts.tearpx, opts.minweight, opts.sigma, round(opts.batch));
if opts.layermaps,   args = [args ' --layermaps'];     end
if ~opts.boundaries, args = [args ' --no-boundaries']; end

cmd = sprintf('"%s" -u %s', opts.pythonexe, args);
if ispc
    % conda DLLs (VTK etc.) resolve only with the env dirs on PATH
    envroot = fileparts(opts.pythonexe);
    pathadd = strjoin({envroot, fullfile(envroot, 'Library', 'bin'), ...
        fullfile(envroot, 'Library', 'mingw-w64', 'bin'), ...
        fullfile(envroot, 'Library', 'usr', 'bin'), ...
        fullfile(envroot, 'Scripts'), fullfile(envroot, 'bin')}, ';');
    cmd = sprintf('set "PATH=%s;%%PATH%%" && %s', pathadd, cmd);
end

if opts.verbose, fprintf('cellsToFlatmap: running Python worker...\n'); end
[status, cmdout] = system(cmd);
if opts.verbose && ~isempty(strtrim(cmdout)), fprintf('%s\n', strtrim(cmdout)); end
if status ~= 0
    error('cellsToFlatmap:pythonFailed', 'Python worker failed (status %d):\n%s', status, cmdout);
end
if exist(outfile, 'file') ~= 2
    error('cellsToFlatmap:noOutput', 'Python worker produced no output. Log:\n%s', cmdout);
end

% ---- collect results (python 0-based pixel -> MATLAB 1-based) -----------
R   = load(outfile);
out = struct();
out.flatxy        = double(R.flatxy) + 1;
out.depth_um      = double(R.depth_um(:));
out.depthnorm_um  = double(R.depthnorm_um(:));
out.layer         = double(R.layer(:));
out.layeredges_um = double(R.layeredges_um);
edges             = double(R.ref_layeredges_um(:)).';
names             = R.layer_names;
if ~iscell(names), names = num2cell(names); end
out.layers = struct('names', {cellfun(@(s) char(string(s)), names(:).', 'UniformOutput', false)}, ...
    'top_um', edges(1:end-1), 'bottom_um', edges(2:end));
out.mapped        = logical(R.mapped(:));
out.isright       = logical(R.isright(:));
out.tear          = logical(R.tear(:));
out.fallback      = logical(R.fallback(:));
out.viewmask      = logical(R.viewmask);
out.midline_x     = double(R.midline_x) + 1;
if isfield(R, 'layermaps_um'), out.layermaps_um = R.layermaps_um; end
out.regionBoundaries = struct('names', {{}}, 'coords', {{}}, 'isright', false(1, 0));
if isfield(R, 'region_names') && isfield(R, 'region_coords')
    nm = R.region_names; co = R.region_coords;
    if ~iscell(nm), nm = num2cell(nm); end
    if ~iscell(co), co = num2cell(co); end
    out.regionBoundaries.names  = cellfun(@(s) char(string(s)), nm(:).', 'UniformOutput', false);
    out.regionBoundaries.coords = cellfun(@(c) double(c) + 1, co(:).', 'UniformOutput', false);
    out.regionBoundaries.isright = logical(R.region_isright(:).');
end
out.coordtype  = opts.coordtype;
if opts.verbose
    fprintf('cellsToFlatmap: %d / %d cells mapped onto the flatmap.\n', ...
        nnz(out.mapped), numel(out.mapped));
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
