function [densout, indsused] = standardAreaDensities()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
dp               = fullfile(allen_atlas_path, 'ero_et_al_2018_densities.CSV');
dpallen          = fullfile(allen_atlas_path, 'parcellation_to_parcellation_term_membership.csv');
densfile         = readtable(dp, 'VariableNamingRule','preserve');
parcelfile       = readtable(dpallen);

ikeep            = strcmp(parcelfile.parcellation_term_set_name,'substructure');
substrtable      = parcelfile(ikeep, :);
substrtable      = sortrows(substrtable, "parcellation_index", 'ascend');
indsused         = substrtable.parcellation_index;
densout         = nan(size(substrtable, 1), 3);
for ii = 1:size(densout, 1)

    nameallen = substrtable.parcellation_term_name{ii};
    nameallen = strrep(lower(nameallen),',','');
    indsfind  = find(contains(lower(densfile.Regions), nameallen));
    
    if numel(indsfind) == 1
        currdens = [densfile.("Neurons [mm-3]")(indsfind) ...
            densfile.("Excitatory [mm-3]")(indsfind) ...
            densfile.("Glia [mm-3]")(indsfind)];
        densout(ii, :) = currdens;
    end
end


end