function finalpts = slicePointsToAtlas(inputpts, trstruct, savepath)
%SLICEPOINTSTOATLAS Map points detected on sections into atlas voxels.
%
%   FINALPTS = SLICEPOINTSTOATLAS(INPUTPTS, TRSTRUCT) takes [x y section]
%   points in pixels of the aligned sections (volume_aligned, px_process; the
%   'cell_locations' convention of extractCellsFromSliceVolume) and returns
%   10 um atlas voxels [DV AP ML] (= [dim2 dim1 dim3] of annotation_10, the
%   order brain cells use). Columns 4:end are carried over.
%
%   The chain is the exact reverse of the one generateRegisteredSliceVolume
%   uses to build volume_registered (mirrored pixel for pixel in
%   sliceAtlasPointsToSection), so a cell lands on the atlas voxel that shows
%   it in the registered volume:
%     section pixel u -> physical (u - 1) * px_process  (elastix origin 0)
%     -> per-section FORWARD B-spline (elastix_forward/NNN_slice_*), which
%        moves section points onto the affinely warped atlas plane
%     -> 10 um plane pixel -> inverse of the atlas-to-image affine
%     -> sample world at the section's AP position -> inverse rigid
%     -> atlas voxel.
%   The reverse B-spline files are for images (see READBSPLINEFIELD).
%
%   FINALPTS = SLICEPOINTSTOATLAS(INPUTPTS, TRSTRUCT, SAVEPATH) reads
%   sliceinfo.mat (pixel sizes) and elastix_forward from SAVEPATH, for
%   registration folders that have moved since transform_params was written.
%   Without it, pixel sizes come from TRSTRUCT (processres, registres,
%   atlasres).
%
%   See also SLICEATLASPOINTSTOSECTION, GENERATEREGISTEREDSLICEVOLUME,
%   READBSPLINEFIELD, APPLYBSPLINEFIELD.

if nargin < 3; savepath = ''; end
px = slicePixelSizes(trstruct, savepath);

%--------------------------------------------------------------------------
% the reference grids sliceAtlasPointsToSection / generateRegisteredSliceVolume use
Ra0 = trstruct.space3d_atlas_20um;
Ra  = imref3d(Ra0.ImageSize * 2, Ra0.XWorldLimits, Ra0.YWorldLimits, Ra0.ZWorldLimits);
Rs0 = trstruct.space3d_sample_20um;
sz2 = Ra.ImageSize(2:3);                                   % in-plane [DV ML] at 10 um
Rs  = imref3d([ceil(Rs0.ImageExtentInWorldY), sz2], ...
    Rs0.XWorldLimits, Rs0.YWorldLimits, Rs0.ZWorldLimits);
sc    = px.register / px.atlas;                            % 10 um px per 20 um unit
regmm = px.register * 1e-3;
sids  = trstruct.sliceids_in_sample_space(:);
Nsec  = numel(sids);
fwdfiles = forwardSliceFiles(trstruct.tformbspline_samp20um_to_atlas_20um, savepath, Nsec);

%--------------------------------------------------------------------------
N      = size(inputpts, 1);
atlas  = nan(N, 3);
secidx = double(inputpts(:, 3));
for k = unique(secidx(secidx >= 1 & secidx <= Nsec))'
    ii = secidx == k;
    m  = double(inputpts(ii, 1:2));                        % [col row] section px
    pm = (m - 1) * px.process * 1e-3;                      % mm, origin 0
    % forward b-spline onto the affinely warped atlas plane (20 um grid file)
    D  = readBsplineField(fwdfiles{k}, regmm);
    q  = (applyBsplineField(D, pm / regmm + 1) - 1) * regmm;   % mm
    b  = q / (px.atlas * 1e-3) + 1;                        % [col row] on the 10 um plane
    % plane -> atlas-plane pixel before the per-section affine
    w  = trstruct.tformaffine_tform_atlas_to_image(k).transformPointsInverse(b / sc);
    pml = w(:, 1) * sc;  cdv = w(:, 2) * sc;
    % section plane -> sample world -> atlas world -> atlas voxel
    [xs, ys, zs] = intrinsicToWorld(Rs, cdv, sids(k) * ones(size(cdv)), pml);
    aw = trstruct.tformrigid_allen_to_samp_20um.transformPointsInverse([xs ys zs]);
    [dv, apr, ml] = worldToIntrinsic(Ra, aw(:, 1), aw(:, 2), aw(:, 3));
    atlas(ii, :) = [dv, apr + trstruct.atlasaplims(1) - 1, ml];
end
%--------------------------------------------------------------------------
finalpts = [atlas, double(inputpts(:, 4:end))];
end

%--------------------------------------------------------------------------
function px = slicePixelSizes(trstruct, savepath)
% pixel sizes (um) of the sections, the registration grid and the atlas
si = struct();
f  = fullfile(char(savepath), 'sliceinfo.mat');
if ~isempty(savepath) && exist(f, 'file')
    s = load(f); if isfield(s, 'sliceinfo'), si = s.sliceinfo; else, si = s; end
end
px.process  = firstField({si, trstruct}, {'px_process', 'processres'}, []);
px.register = firstField({si, trstruct}, {'px_register', 'registres'}, 20);
px.atlas    = firstField({si, trstruct}, {'px_atlas', 'atlasres'}, 10);
if isempty(px.process)
    error('slicePointsToAtlas:noPixelSize', ['Section pixel size unknown: pass the ' ...
        'registration folder (it holds sliceinfo.mat) as the third input.']);
end
end

function v = firstField(structs, names, default)
v = default;
for s = 1:numel(structs)
    for n = 1:numel(names)
        if isfield(structs{s}, names{n}) && ~isempty(structs{s}.(names{n}))
            v = double(structs{s}.(names{n})); return;
        end
    end
end
end

%--------------------------------------------------------------------------
function files = forwardSliceFiles(revdir, savepath, Nslices)
% the forward files sit in elastix_forward: in savepath if given, else next
% to the stored reverse folder
revdir  = regexprep(char(revdir), '[\\/]+$', '');
fwddirs = {fullfile(fileparts(revdir), 'elastix_forward')};
if ~isempty(savepath)
    fwddirs = [{fullfile(char(savepath), 'elastix_forward')}, fwddirs];
end
ifound = find(cellfun(@(d) exist(d, 'dir') == 7, fwddirs), 1);
if isempty(ifound)
    error('slicePointsToAtlas:noForwardFile', ...
        ['No elastix_forward folder found (looked in: %s). Points need the ' ...
         'forward transforms that registerSlicesToAtlas writes there; pass ' ...
         'the registration folder as the third input if it has moved.'], ...
        strjoin(fwddirs, ', '));
end
fwddir = fwddirs{ifound};
files  = cell(Nslices, 1);
for islice = 1:Nslices
    files{islice} = fullfile(fwddir, sprintf('%03d_slice_bspline_atlas_to_samp_20um.txt', islice));
    if ~exist(files{islice}, 'file')
        error('slicePointsToAtlas:noForwardFile', 'Forward B-spline file %s not found.', files{islice});
    end
end
end
