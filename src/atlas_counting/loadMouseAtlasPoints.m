function [atlaspts,filetoload] = loadMouseAtlasPoints(mouseid, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    pointssample = false;
else
    pointssample = varargin{1};
end

if pointssample
    filelook  = 'cell_locations_sample.mat';
    namecheck = 'cell_locations';
else
    filelook = 'cell_locations_atlas.mat';
    namecheck = 'atlasptcoords';
end

savepath  = 'D:\DATA_folder\Mice';
toget     = dir(fullfile(savepath, mouseid, 'Anatomy', filelook));
if ~isempty(toget)
    filetoload = fullfile(toget.folder, toget.name);
    atlaspts   = load(filetoload, namecheck);
    atlaspts   = atlaspts.(namecheck);
else
    atlaspts   = [];
    filetoload = '';
end
end