function [tv, av, pinfo, atlasres] = loadAtlasInfo(atlasid)
%LOADATLASINFO Template, annotation and region table of a supported atlas.
%
%   [TV, AV, PINFO, ATLASRES] = LOADATLASINFO(ATLASID) with ATLASID one of
%       'allen2020_10um'      Allen CCFv3 (2020 parcellation), 10 um
%       'allen2020fusi_50um'  vessel atlas for fUS, 50 um
%       'spinalcord'          spinal cord atlas (Fiederling et al. 2021),
%                             native 10 x 10 x 20 um
%   ATLASRES is the voxel size in mm. PINFO is the region table: the Allen
%   parcellation_to_parcellation_term_membership.csv, or Atlas_Regions.csv for
%   the cord. Use atlasNameLookup(PINFO) to get region names from AV values.

loadav = nargout > 1;

switch atlasid
    case 'allen2020_10um'
        allen2020path = fileparts(which('annotation_10.nii.gz'));
        if isempty(allen2020path)
            error('No CCF atlas found (add CCF atlas to path)')
        end
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
    case 'spinalcord'
        % volumes at native resolution, as distributed (no long-axis flip)
        [tv, av, pinfo, ~, atlasresum] = loadSpinalCordAtlas();
        atlasres = atlasresum * 1e-3;
    case 'youratlas'
        % add your own...
        % tv = ...
        % av = ...
        % pinfo = ...
    otherwise
        error('loadAtlasInfo:unknownAtlas', 'Unknown atlas id "%s".', atlasid);
end
end
