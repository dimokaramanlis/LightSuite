function finalpts = slicePointsToAtlas(inputpts, trstruct, savepath)
%SLICEPOINTSTOATLAS Map points detected on sections into atlas voxels.
%
%   FINALPTS = SLICEPOINTSTOATLAS(INPUTPTS, TRSTRUCT) takes [x y slice] points
%   (columns 4:end are carried over) and runs, per section: B-spline onto the
%   affinely warped atlas plane -> inverse of the atlas-to-image affine ->
%   rigid allen_to_samp inverse -> 10 um atlas voxels.
%
%   The B-spline step uses each section's FORWARD file
%   (elastix_forward/NNN_slice_bspline_atlas_to_samp_20um.txt), which moves
%   section points onto the atlas plane, and adds its displacement; the
%   reverse files in TRSTRUCT are for images (see READBSPLINEFIELD).
%
%   FINALPTS = SLICEPOINTSTOATLAS(INPUTPTS, TRSTRUCT, SAVEPATH) looks for
%   SAVEPATH/elastix_forward first, for registration folders that have moved
%   since transform_params.mat was written.
%
%   See also READBSPLINEFIELD, APPLYBSPLINEFIELD, REGISTERSLICESTOATLAS.

%--------------------------------------------------------------------------
Nslices    = numel(trstruct.sliceids_in_sample_space);
regsize_mm = trstruct.atlasres * 2 * 1e-3;
if nargin < 3; savepath = ''; end
fwdfiles   = forwardSliceFiles(trstruct.tformbspline_samp20um_to_atlas_20um, savepath, Nslices);
outpts     = inputpts;
Rsample3d  = trstruct.space3d_sample_20um;
ysampvals  = linspace(Rsample3d.YWorldLimits(1)+0.5, Rsample3d.YWorldLimits(2)-0.5,...
    Rsample3d.ImageExtentInWorldY);

for islice = 1:Nslices
    %----------------------------------------------------------------------
    currids  = inputpts(:, 3) == islice;
    if ~any(currids); continue; end
    currlocs = double(inputpts(currids, 1:3));
    pts      = currlocs(:, 1:2) * 1.25 * 1e-3/regsize_mm;   % 20 um voxels [x y]
    afftrans = trstruct.tformaffine_tform_atlas_to_image(islice);
    %----------------------------------------------------------------------
    % forward b-spline onto the affinely warped atlas plane
    Dfield   = readBsplineField(fwdfiles{islice}, regsize_mm);
    finalpts = applyBsplineField(Dfield, pts, 'forward');
    %----------------------------------------------------------------------
    % we then invert affines
    finalpts = afftrans.transformPointsInverse(finalpts);
    outpts(currids, 1:2) = finalpts;
    outpts(currids,   3) = ysampvals(trstruct.sliceids_in_sample_space(islice));
    %----------------------------------------------------------------------
end
%-------------------------------------------------------------------------
rigidfull          = trstruct.tformrigid_allen_to_samp_20um.invert;
atlasptcoords      = rigidfull.transformPointsForward(outpts(:, [2 3 1]))*2;
atlasptcoords(:,2) = atlasptcoords(:,2) + trstruct.atlasaplims(1);
%-------------------------------------------------------------------------
finalpts = [atlasptcoords inputpts(:, 4:end)];
%-------------------------------------------------------------------------
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
