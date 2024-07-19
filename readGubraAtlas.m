function [tv,av, gubra_atlas_path] = readGubraAtlas()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gubra_atlas_path = fileparts(which('gubra_template_olf.nii.gz'));
if isempty(gubra_atlas_path)
    error('No atlas found (add atlas to path)')
end

disp('Loading Gubra atlas...')
volsread = {'gubra_template_olf.nii.gz', 'gubra_ano_olf.nii.gz'};

tv = prepareVolume(fullfile(gubra_atlas_path, volsread{1}));
av = prepareVolume(fullfile(gubra_atlas_path, volsread{2}));


end

function Vatlas = prepareVolume(dp)
Vatlas  = niftiread(dp);
Vatlas  = permute(Vatlas, [2 3 1]);
Vatlas  = flip(Vatlas, 3);
Vatlas  = flip(Vatlas, 2);
end