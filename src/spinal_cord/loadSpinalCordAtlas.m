function [tv, av, parcelinfo, segmentinfo] = loadSpinalCordAtlas()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%-------------------------------------------------------------------------
dplook   = fileparts(which("Segments.csv"));
dpav     = fullfile(dplook, "Annotation.tif");
dptv     = fullfile(dplook, "Template.tif");
tv       = tiffreadVolume(dptv);
av           = tiffreadVolume(dpav);
atlasres     = [10 10 20];
parcelinfo   = readtable(fullfile(dplook,'Atlas_Regions.csv'));
segmentinfo  = readtable(fullfile(dplook,'Segments.csv'));
%-------------------------------------------------------------------------
end