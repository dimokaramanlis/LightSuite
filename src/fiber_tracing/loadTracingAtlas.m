function atlas = loadTracingAtlas()
% LOADTRACINGATLAS  Load the Allen CCF annotation volume and parcellation
%   look-up used by the tracing tools (GRIN lens / Neuropixels probes).
%
%   ATLAS = loadTracingAtlas()
%
%   Requires annotation_10.nii.gz and
%   parcellation_to_parcellation_term_membership.csv to be on the MATLAB path
%   (the standard LightSuite atlas resources).
%
%   Output struct ATLAS
%   -------------------
%   .av        – 3-D annotation volume of parcellation indices, indexed as
%                av(AP, DV, ML).  Outside-brain voxels are 0.
%   .areaidx   – column vector of unique substructure parcellation indices.
%   .acronyms  – cell array of region acronyms, aligned with .areaidx.
%   .names     – cell array of full region names, aligned with .areaidx.
%   .colors    – N×3 RGB (0–255) region colours, aligned with .areaidx.

    allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
    if isempty(allen_atlas_path)
        error('loadTracingAtlas: annotation_10.nii.gz not found on the MATLAB path.');
    end
    atlas.av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));

    parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));
    if isempty(parcel_path)
        error('loadTracingAtlas: parcellation_to_parcellation_term_membership.csv not found on the MATLAB path.');
    end
    parcelinfo = readtable(fullfile(parcel_path, ...
        'parcellation_to_parcellation_term_membership.csv'));

    % Restrict to the 'substructure' term set – the level stored in annotation_10
    substridx     = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib] = unique(parcelinfo.parcellation_index(substridx));

    acronyms = parcelinfo.parcellation_term_acronym(substridx);
    names    = parcelinfo.parcellation_term_name(substridx);
    colors   = [parcelinfo.red(substridx), parcelinfo.green(substridx), parcelinfo.blue(substridx)];

    atlas.areaidx  = areaidx;
    atlas.acronyms = acronyms(ib);
    atlas.names    = names(ib);
    atlas.colors   = colors(ib, :);
end
