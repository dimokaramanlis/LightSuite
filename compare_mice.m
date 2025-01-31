mouseIds = {'DK025', 'DK027', 'DK026', 'DK028', 'DK031'};



% DK025, DK027 joint

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
%%

Ngroups = max(av, [], 'all');
Nmice   = numel(mouseIds);
groupstrs       = lower(table2cell(st(1:Ngroups,4)));

countsall = nan(Ngroups, Nmice);

for mouseid = 1:Nmice
    atlasptcoords = loadMouseAtlasPoints(mouseIds{mouseid});
    cleanatlaspts                    = sanitizeCellCoords(atlasptcoords, av);
    [areacounts, areavols, idcareas] = groupCellsIntoLeafRegions(cleanatlaspts, av);
    countsall(:, mouseid) = areacounts;
%     avplot = ndSparse.build(cleanatlaspts(:,[2 1 3]), 1,size(av));
% avplot = single(full(avplot));
% avgauss = imgaussfilt3(avplot, 10);

end

%%

irem = sum(isnan(densitiesall),2)==Nmice;
groupstrsuse = groupstrs(~irem, :);
densitiesuse = countsall(~irem, :);
volsuse      = areavols(~irem, :);

summotor     = sum(densitiesuse(contains(groupstrsuse,'primary motor'), :), 1);
densitiesuse = (densitiesuse./volsuse)./summotor;

ivis = contains(groupstrsuse,'primary motor')|...
    contains(groupstrsuse,'primary visual')|...
    contains(groupstrsuse,'prelimbic')|...
    contains(groupstrsuse,'infralimbic')|...
    contains(groupstrsuse,'cingulate');

% ivis = contains(groupstrsuse,'cingulate')|...
%     contains(groupstrsuse,'primary motor');

bar(groupstrsuse(ivis), densitiesuse(ivis,:))
ylabel('Density (norm.)')
legend(mouseIds)%%
%%
socidx     = (mean(densitiesuse(:, 1:2), 2) - mean(densitiesuse(:, 3:4), 2))./(mean(densitiesuse(:, 1:2), 2) + mean(densitiesuse(:, 3:4), 2));
[~, isort] = sort(socidx, 'descend');

plot(densitiesuse(isort,:))




%%
%%
stmat                            = atlasTreeToMat(st);
Nmax = numel(areacounts);
%%
areavolsuse     = areavols; % in mm3
cells_per_group = areacounts;
irem = cells_per_group<2;
cells_per_group(irem) = [];
areavolsuse(irem) = [];
stgroups = stmat(1:Nmax, :);

% % move one up
% lastleaf   = findfirst(~isnan(stgroups),2,1,'last');
% prevbranch = lastleaf - 1;
% idxprev    = sub2ind(size(stgroups), find(prevbranch>0), prevbranch(prevbranch>0));
% 
% 
% idxprev - 1

groupstrs(irem) = [];
