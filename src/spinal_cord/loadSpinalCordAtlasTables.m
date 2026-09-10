function [parcelinfo, segmentinfo, atlasres] = loadSpinalCordAtlasTables()
%LOADSPINALCORDATLASTABLES Region and segment tables of the spinal cord atlas.
%
%   [PARCELINFO, SEGMENTINFO] = LOADSPINALCORDATLASTABLES() reads
%   Atlas_Regions.csv and Segments.csv without touching the template and
%   annotation volumes. Analysis scripts that work from the per-region
%   statistics LightSuite already saved need the names and the hierarchy, not
%   the volumes, and those are a gigabyte of reading.
%
%   [..., ATLASRES] = LOADSPINALCORDATLASTABLES() also returns the native voxel
%   size, [10 10 20] um. Segment start and end planes are indices into the
%   native grid.
%
%   See also LOADSPINALCORDATLAS, REORGANIZESPINALCORDAREAS.

dplook = fileparts(which('Segments.csv'));
if isempty(dplook)
    error('loadSpinalCordAtlasTables:atlasNotFound', ...
        ['Could not find the spinal cord atlas: Segments.csv is not on the ' ...
         'MATLAB path. Add the atlas folder to the path and try again.']);
end

parcelinfo  = readtable(fullfile(dplook, 'Atlas_Regions.csv'));
segmentinfo = readtable(fullfile(dplook, 'Segments.csv'));
atlasres    = [10 10 20];

end
