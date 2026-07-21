function [tv, av, pinfo, atlasres] = loadAtlasInfo(atlasid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

loadav = nargout > 1;

switch atlasid
    case 'allen2020_10um'
        allen2020path = fileparts(which('annotation_10.nii.gz'));
        pinfo         = readtable(fullfile(allen2020path, 'parcellation_to_parcellation_term_membership.csv'));
        tv            = niftiread(fullfile(allen2020path,"average_template_10.nii.gz"));
        atlasres      = [1 1 1]*0.01;
        if loadav
            av  = niftiread(fullfile(allen2020path,"annotation_10.nii.gz"));
        end
    case 'allen2020fusi_50um'
        allen2020path  = fileparts(which('annotation_10.nii.gz'));
        fusiatlaspath  = 'D:\fusi_test\new_atlas';
        tv       = niftiread(fullfile(fusiatlaspath, 'vessel_atlas_template.nii.gz'));
        if loadav
            av       = niftiread(fullfile(fusiatlaspath,'vessel_atlas_annotation.nii.gz'));
        end
        pinfo     = readtable(fullfile(allen2020path, 'parcellation_to_parcellation_term_membership.csv'));
        atlasres      = [1 1 1]*0.05;
    case 'youratlas'
        % add your own...
        % tv = ...
        % av = ...
        % pinfo = ...
end
end