function [D, tformpath] = readBsplineField(source, spacing_mm, savepath)
%READBSPLINEFIELD Displacement field of an elastix B-spline file, in voxels.
%
%   D = READBSPLINEFIELD(TFORMPATH, SPACING_MM) runs transformix on the
%   transform file and returns its displacement field divided by the voxel
%   spacing: [Sx Sy Sz 3] in 3D, [Sx Sy 2] in 2D, with x = image columns,
%   y = rows, z = planes. D(i,j,k,:) is the displacement of the 1-based voxel
%   (i,j,k), whose physical position is ((i,j,k) - 1) * SPACING_MM, because
%   LightSuite always hands elastix images with origin 0. Move points with
%   APPLYBSPLINEFIELD.
%
%   D = READBSPLINEFIELD(TRSTRUCT, SPACING_MM, SAVEPATH) takes the
%   transform_params struct of a brain or spinal cord registration and reads
%   its FORWARD file, bspline_atlas_to_samp_20um.txt, the file that moves
%   sample points onto the atlas. transform_params only stores the path of
%   the inverse file (bspline_samp_to_atlas_*), which is for images; the
%   forward file is saved in the same folder. SAVEPATH (optional) is searched
%   first, so a registration folder that has been moved still works.
%
%   Which file for which job (elastix moves points and images in opposite
%   directions):
%       sample points -> atlas     forward file, field added
%       atlas images  -> sample    forward file
%       sample images -> atlas     inverse file
%       atlas points  -> sample    forward file, inverted (APPLYBSPLINEFIELD)
%
%   [D, TFORMPATH] = ... also returns the file that was read.
%
%   See also APPLYBSPLINEFIELD, BRAINPOINTTRANSFORM, CORDPOINTSTOATLAS.

if nargin < 3
    savepath = '';
end

if isstruct(source)
    tformpath = forwardFileFromParams(source, savepath);
else
    tformpath = char(source);
end
if ~exist(tformpath, 'file')
    error('readBsplineField:missingFile', 'B-spline file not found: %s', tformpath);
end

F    = transformix([], tformpath, 'verbose', false);   % [ndim x Sx Sy (Sz)], mm
ndim = size(F, 1);
D    = permute(F, [2:ndim+1, 1]) / spacing_mm;
end

%--------------------------------------------------------------------------
function fwdpath = forwardFileFromParams(trstruct, savepath)
fwdname    = 'bspline_atlas_to_samp_20um.txt';
candidates = {};
if ~isempty(savepath)
    candidates{end+1} = fullfile(char(savepath), fwdname);
end
if isfield(trstruct, 'tform_bspline_samp20um_to_atlas_20um_px')
    candidates{end+1} = fullfile(fileparts(trstruct.tform_bspline_samp20um_to_atlas_20um_px), fwdname);
end
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file')
        fwdpath = candidates{k};
        return;
    end
end
error('readBsplineField:noForwardFile', ...
    ['Forward B-spline file %s not found (looked in: %s). Points need the ' ...
     'forward transform; it is written by the B-spline registration next to ' ...
     'the inverse file.'], fwdname, strjoin(candidates, ', '));
end
