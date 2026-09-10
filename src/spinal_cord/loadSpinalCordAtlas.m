function [tv, av, parcelinfo, segmentinfo, atlasres, nativesize] = loadSpinalCordAtlas(outputres)
%LOADSPINALCORDATLAS Load the spinal cord atlas template, annotation and tables.
%
%   [TV, AV, PARCELINFO, SEGMENTINFO] = LOADSPINALCORDATLAS() returns the
%   template volume TV and the annotation volume AV at the native atlas
%   resolution, together with the region table PARCELINFO (Atlas_Regions.csv)
%   and the segment table SEGMENTINFO (Segments.csv). The atlas is found on the
%   MATLAB path through Segments.csv.
%
%   [...] = LOADSPINALCORDATLAS(OUTPUTRES) resamples both volumes to the voxel
%   size OUTPUTRES (1x3, in micrometres, per array dimension [row col plane]),
%   which is how the registration pipeline pulls the atlas in at its own
%   resolution. The annotation is resampled with nearest-neighbour so label
%   values stay intact.
%
%   [..., ATLASRES, NATIVESIZE] = LOADSPINALCORDATLAS(...) also returns the
%   *native* voxel size of the atlas, [10 10 20] um, and the size of the native
%   volumes, whatever OUTPUTRES was. Segment start and end planes in SEGMENTINFO
%   are indices into that native grid.
%
%   See also LOADSPINALCORDATLASANDPOINTS, LOADCORDATLASVOLUMES.

%--------------------------------------------------------------------------
dplook = fileparts(which('Segments.csv'));
if isempty(dplook)
    error('loadSpinalCordAtlas:atlasNotFound', ...
        ['Could not find the spinal cord atlas: Segments.csv is not on the ' ...
         'MATLAB path. Add the atlas folder to the path and try again.']);
end
%--------------------------------------------------------------------------
tv          = tiffreadVolume(fullfile(dplook, 'Template.tif'));
av          = tiffreadVolume(fullfile(dplook, 'Annotation.tif'));
atlasres    = [10 10 20];
parcelinfo  = readtable(fullfile(dplook, 'Atlas_Regions.csv'));
segmentinfo = readtable(fullfile(dplook, 'Segments.csv'));
nativesize  = size(av);
%--------------------------------------------------------------------------
if nargin > 0 && ~isempty(outputres)
    scalefac = atlasres ./ reshape(outputres, 1, 3);
    if ~all(scalefac == 1)
        tv = imresize3(tv, 'Scale', scalefac);
        av = imresize3(av, 'nearest', 'Scale', scalefac);
    end
end
%--------------------------------------------------------------------------
end
