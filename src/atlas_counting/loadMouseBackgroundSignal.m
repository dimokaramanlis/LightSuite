function backsig = loadMouseBackgroundSignal(mouseid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

savepath  = 'D:\DATA_folder\Mice';
toget     = dir(fullfile(savepath, mouseid, 'Anatomy', 'background_volume_areas.mat'));
backsig   = load(fullfile(toget.folder, toget.name));
if isfield(backsig, 'backvolareas')
    backsig   = backsig.backvolareas;
else
    backsig   = backsig.medianoverareas;
end
backsig(backsig<0) = nan;
end