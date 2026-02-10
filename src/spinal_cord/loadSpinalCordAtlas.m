function [tv, av, parcelinfo, segments] = loadSpinalCordAtlas(dplook, outputres)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%-------------------------------------------------------------------------
dpav     = fullfile(dplook, "annotation.tiff");
dptv     = fullfile(dplook, "reference.tiff");
tv       = tiffreadVolume(dptv);
av       = tiffreadVolume(dpav);
atlasres = [10 10 20];
%-------------------------------------------------------------------------
parcelinfo  = readtable(fullfile(dplook,'structures.csv'));

segments    = readtable(fullfile(dplook,'Segments.csv'));

%-------------------------------------------------------------------------

tv       = imresize3(tv, 'Scale', atlasres./outputres);
av       = imresize3(av, 'nearest', 'Scale', atlasres./outputres);

%-------------------------------------------------------------------------
end