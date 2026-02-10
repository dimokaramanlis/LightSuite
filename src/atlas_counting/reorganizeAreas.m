function resout = reorganizeAreas(counts, signals, volumes, parcelinfo, aggtype)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% count and volume correspond to half a hemisphere
%--------------------------------------------------------------------------
[Ngroups, ~, Nmice] = size(counts);
% smat    = atlasTreeToMat(st);
% smat    = smat(1:Ngroups, :);
[agglist, aggnames, parcelout] = getAggList(parcelinfo, aggtype);
[agglistcoarse, coarsenames]   = getAggList(parcelinfo, 'division');

indsout = parcelout.parcellation_index;
%--------------------------------------------------------------------------
% we first organize areas in a reasonable depth

cout       = nan(numel(aggnames), Nmice);
sigout     = nan(numel(aggnames), Nmice);
volout     = nan(numel(aggnames), 1);
nameout    = cell(numel(aggnames), 3);
for ii = 1:numel(aggnames)
    row  = find(agglist== ii);
    % counts are averaged in hemispheres
    countsignal = squeeze(mean(counts(row, :, :), 2, 'omitmissing'));
    % countsignal = squeeze(max(counts(row, :, :), [], 2, 'omitmissing'));

    countsignal = reshape(countsignal, numel(row), Nmice);
    cout(ii, :) = sum(countsignal, 1, 'omitmissing');
    % volumes are the whole volme
    volssel         = volumes(row);
    volout(ii, :)   = sum(volssel, 1, 'omitmissing')/2;
    
    % signal is averaged based on volume area
    wtscurr      = volssel/sum(volssel);
    bkgsignal    = squeeze(mean(signals(row, :, :), 2, 'omitmissing'));
    bkgsignal    = reshape(bkgsignal, numel(row), Nmice);
    signaluse    = sum(wtscurr.*bkgsignal, 1, 'omitmissing');
    sigout(ii,:) = signaluse;


    lastname = coarsenames(agglistcoarse(row));
    nameout(ii, :) = {parcelout.parcellation_term_name{ii}, ...
        parcelout.parcellation_term_acronym{ii}, lastname{1}};

end

%--------------------------------------------------------------------------
% we then get rid of things 
% 
% 1. without volume in the atlas
% 2. Hindbrain, Cerebellum, Olfactory bulb, fiber tracts

iclassrem = contains(lower(nameout(:,3)), 'nerves') |...
    contains(lower(nameout(:,3)), 'tracts') |...
    contains(lower(nameout(:,3)), 'system') | ...
    contains(lower(nameout(:,3)), 'unassigned') | ...
    contains(lower(nameout(:,3)), 'white matter') | ...
    contains(lower(nameout(:,3)), 'canal') | ...
    contains(lower(nameout(:,3)), 'aqueduct') | ...
    contains(lower(nameout(:,3)), 'ventricle');

% ikeeptype   = ~any(ismember(smat, st.id(find(iclassrem))), 2);
% ikeepvolume = volout > 1e-5;

ikeep   = ~iclassrem; % & ikeeptype;
resout.counts   = cout(ikeep,  :);
resout.volumes  = volout(ikeep, :);
resout.signal   = sigout(ikeep, :);
resout.names    = nameout(ikeep, :);
resout.indices  = indsout(ikeep, :);
% smatnew    = smat(ikeep, :);

%--------------------------------------------------------------------------
% areaids   = sorg(toindex);
% [unids, ~, ic] = unique(areaids);
% Nareas         = numel(unids);
% 
% finalvolumes = accumarray(ic, volumesnew, [Nareas 1], @sum);
% finalcounts  = nan(Nareas, 2, Nmice, 'single');
% for imouse = 1:Nmice
%     for ii = 1:2 
%         finalcounts(:, ii, imouse) = accumarray(ic, countsnew(:, ii, imouse), [Nareas 1], @nansum);
%     end
% end
% finalids = lower(st.name(findfirst(unids == st.id', 2)));
% %--------------------------------------------------------------------------
% 
% ivis = contains(finalids,'primary motor')|...
%     contains(finalids,'primary visual')|...
%     contains(finalids,'prelimbic')|...
%     contains(finalids,'infralimbic')|...
%     contains(finalids,'primary auditory')|...
%     contains(finalids,'secondary motor')|...
%     contains(finalids,'cingulate');
% 
% 
% tomerge  = find(~isnan(smat(:, depth + 1)));
% lowdepth = find(isnan(smat(:, depth + 1)));
% %--------------------------------------------------------------------------
% % we first keep the lowdepth areas
% countslowdepth  = counts(lowdepth(volumes(lowdepth) > 0), :, :);
% volumeslowdepth = volumes(lowdepth(volumes(lowdepth) > 0));
% 
% 
% %--------------------------------------------------------------------------
% 
% smat()
% 
% [newcheck,~,ic] = unique(smat(:,depth));
% ikeep   = ~isnan(newcheck);
% newvols = accumarray(ic, volumes, [], @sum);
% newvols = newvols(ikeep);
% 
% locations = findIndicesLocation(newcheck(ikeep), st.id(1:Ngroups));
% 
% st.name(locations)
% 
% ismember(st.id(1:Ngroups), newcheck(ikeep))
% 

end
    
function [aggList, aggNames, parcelout] = getAggList(parcelinfo, aggType)

ichoose       = strcmp(parcelinfo.parcellation_term_set_name, aggType);
reducedparcel = parcelinfo(ichoose, :);
reducedparcel = sortrows(reducedparcel, "parcellation_index", 'ascend');

[aggNames, ia, aggList] = unique(reducedparcel.parcellation_term_name);
parcelout = reducedparcel(ia,:);
end