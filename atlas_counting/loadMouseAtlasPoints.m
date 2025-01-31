function atlaspts = loadMouseAtlasPoints(mouseid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

savepath  = 'D:\DATA_folder\Mice';
toget     = dir(fullfile(savepath, mouseid, '**', 'cell_locations_atlas.mat'));
atlaspts  = load(fullfile(toget.folder, toget.name));
atlaspts  = atlaspts.atlasptcoords;
end