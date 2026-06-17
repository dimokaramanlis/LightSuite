function [distFromMidline] = getAllenDistFromMidline()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
parcelinfo       = readtable(fullfile(allen_atlas_path, 'parcellation_to_parcellation_term_membership.csv'));
areaidx          = unique(parcelinfo.parcellation_index);
Ngroups          = numel(areaidx);
Nforaccum        = max(av, [], 'all') + 1;

Npxlr          = size(av,3)/2;
[~,~,zz]       = meshgrid(1:size(av,1), 1:size(av,2),1:Npxlr);
sideav         = reshape(av(:, :, 1:Npxlr), [], 1);
ikeep     = sideav>0;
distFromMidline = accumarray(sideav(ikeep)+1, zz(ikeep), [Nforaccum 1], @median);
distFromMidline = distFromMidline(areaidx+1);
distFromMidline = (Npxlr - distFromMidline)*0.01;


end