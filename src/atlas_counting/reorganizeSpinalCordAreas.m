function resout = reorganizeSpinalCordAreas(counts, signals, volumes, parcelinfo, avinds, aggtype)
%REORGANIZESPINALCORDAREAS Aggregate cord data to a coarser anatomical level.
%
%   RESOUT = REORGANIZESPINALCORDAREAS(COUNTS, SIGNALS, VOLUMES, PARCELINFO,
%   AVINDS, AGGTYPE) rolls per-region cord statistics up to the level AGGTYPE,
%   which is one of
%
%     'substructure' the labelled regions themselves, renamed but not merged
%     'structure'    laminae I-X, and the df, lf and vf funiculi
%     'division'     gray matter and white matter
%
%   the same three levels the Allen brain parcellation uses, so cord results can
%   be read next to brain ones (see REORGANIZEAREAS). CORDAREAHIERARCHY does the
%   anatomy; this function does the arithmetic.
%
%   INPUTS
%     counts     Nareas x Nsegments x Nmice cell counts, or [] when only
%                intensities are being aggregated.
%     signals    Nareas x Nsegments x Nmice intensities, or [] when only counts
%                are being aggregated.
%     volumes    Nareas x Nsegments x Nmice region volumes in mm^3. A
%                Nareas x Nsegments array is taken to hold for every mouse,
%                which is what the atlas volumes are.
%     parcelinfo Atlas_Regions.csv as a table, from LOADSPINALCORDATLASTABLES.
%     avinds     Nareas x 1 atlas ids labelling the rows of the arrays above.
%     aggtype    the level to aggregate to.
%
%   Every level is aggregated the way the quantity itself demands: counts and
%   volumes are summed over the regions that make up a target, while signals are
%   averaged over them weighted by volume, so a large region is not outvoted by
%   a sliver of a neighbour. The rostrocaudal axis is left alone throughout -
%   segments stay separate columns, because a cord region runs the whole length
%   of the cord and one number for it would hide the axis that matters. The
%   segment count never enters the arithmetic, so an atlas with 34 segments and
%   one with any other number both work unchanged.
%
%   OUTPUT
%     RESOUT with fields
%       counts   Ntargets x Nsegments x Nmice, summed  (NaN if COUNTS was [])
%       signal   Ntargets x Nsegments x Nmice, volume-weighted mean
%       volumes  Ntargets x Nsegments x Nmice, summed
%       names    Ntargets x 3 cellstr: name, acronym, division it belongs to
%       indices  Ntargets x 1 atlas ids of the targets
%       level    the AGGTYPE that produced this
%     Targets come in anatomical order - laminae I-X, then df, lf, vf; GM, then
%     WM - and a target no region in AVINDS reaches is dropped.
%
%   Example
%     [parcelinfo, segments] = loadSpinalCordAtlasTables();
%     s = load('chan01_cellcounts.mat');   % areacounts, areavols, areaidx
%     res = reorganizeSpinalCordAreas(s.areacounts, [], s.areavols, ...
%         parcelinfo, s.areaidx, 'structure');
%
%   See also CORDAREAHIERARCHY, REORGANIZEAREAS, LOADSPINALCORDATLASTABLES,
%   CORDATLASGROUPING.

%--------------------------------------------------------------------------
% 1. sizes, from whichever quantity was actually supplied
%--------------------------------------------------------------------------
Nareas = numel(avinds);

reference = counts;
if isempty(reference); reference = signals; end
if isempty(reference); reference = volumes; end
if isempty(reference)
    error('reorganizeSpinalCordAreas:noData', ...
        'At least one of counts, signals or volumes must be non-empty.');
end

Nsegs = size(reference, 2);
Nmice = size(reference, 3);

hascounts  = ~isempty(counts);
hassignals = ~isempty(signals);

counts  = expandToMice(counts,  Nareas, Nsegs, Nmice, 'counts');
signals = expandToMice(signals, Nareas, Nsegs, Nmice, 'signals');
volumes = expandToMice(volumes, Nareas, Nsegs, Nmice, 'volumes');

%--------------------------------------------------------------------------
% 2. which target every row belongs to
%--------------------------------------------------------------------------
[hier, levels] = cordAreaHierarchy(parcelinfo, avinds);

switch lower(char(aggtype))
    case 'substructure'
        % the finest grain is already the answer; only regions that sit inside
        % the hierarchy are kept, so the levels stay nested
        inside  = ~isnan(hier.division_index);
        targets = table(hier.parcellation_index(inside), hier.acronym(inside), ...
            hier.name(inside), 'VariableNames', {'index', 'acronym', 'name'});
        rowlevel = hier.parcellation_index;
    case 'structure'
        targets  = levels.structure;
        rowlevel = hier.structure_index;
    case 'division'
        targets  = levels.division;
        rowlevel = hier.division_index;
    otherwise
        error('reorganizeSpinalCordAreas:badAggType', ...
            'aggtype must be ''substructure'', ''structure'' or ''division'', not ''%s''.', ...
            char(aggtype));
end

%--------------------------------------------------------------------------
% 3. aggregate
%--------------------------------------------------------------------------
Ntargets = height(targets);

cout    = nan(Ntargets, Nsegs, Nmice);
sigout  = nan(Ntargets, Nsegs, Nmice);
volout  = nan(Ntargets, Nsegs, Nmice);
nameout = cell(Ntargets, 3);
indsout = nan(Ntargets, 1);

for ii = 1:Ntargets
    rows = find(rowlevel == targets.index(ii));
    if isempty(rows)
        continue
    end

    volblock = volumes(rows, :, :);
    totalvol = sum(volblock, 1, 'omitnan');

    cout(ii, :, :)   = sum(counts(rows, :, :), 1, 'omitnan');
    volout(ii, :, :) = totalvol;

    % volume-weighted mean, computed per segment and per mouse so that a region
    % missing from one segment does not skew the others
    sigblock = signals(rows, :, :);
    hassig   = ~isnan(sigblock);
    weights  = volblock;
    weights(~hassig | isnan(weights)) = 0;

    % a region with a signal but no volume to weight it by still counts, it just
    % counts equally: without this a caller who has no volumes gets nothing back
    novolume = sum(weights, 1) == 0 & any(hassig, 1);
    weights(repmat(novolume, numel(rows), 1, 1) & hassig) = 1;

    % renormalise over the regions that actually carry a signal here; where none
    % does, the target has no signal in that segment rather than a signal of 0
    wtotal = sum(weights, 1);
    nosig  = wtotal == 0;
    wtotal(nosig) = 1;

    weighted = sum(sigblock .* (weights ./ wtotal), 1, 'omitnan');
    weighted(nosig) = NaN;
    sigout(ii, :, :) = weighted;

    indsout(ii)    = targets.index(ii);
    nameout(ii, :) = {targets.name{ii}, targets.acronym{ii}, ...
        divisionOfTarget(hier, rows)};
end

%--------------------------------------------------------------------------
% 4. drop targets no region reached
%--------------------------------------------------------------------------
keep = ~isnan(indsout);

% a quantity that was never supplied stays missing rather than summing to zero
if ~hascounts;  cout   = nan(size(cout));   end
if ~hassignals; sigout = nan(size(sigout)); end

resout.counts  = cout(keep, :, :);
resout.signal  = sigout(keep, :, :);
resout.volumes = volout(keep, :, :);
resout.names   = nameout(keep, :);
resout.indices = indsout(keep);
resout.level   = lower(char(aggtype));

end

%==========================================================================
% Local helpers
%==========================================================================
function out = expandToMice(in, Nareas, Nsegs, Nmice, label)
%EXPANDTOMICE Bring one quantity to Nareas x Nsegments x Nmice.
%   An empty input becomes all-NaN, so a caller with only intensities or only
%   counts does not have to fabricate the other.

if isempty(in)
    out = nan(Nareas, Nsegs, Nmice);
    return
end

out = double(in);
if size(out, 1) ~= Nareas
    error('reorganizeSpinalCordAreas:sizeMismatch', ...
        '%s has %d rows but avinds lists %d areas.', label, size(out, 1), Nareas);
end
if size(out, 2) ~= Nsegs
    error('reorganizeSpinalCordAreas:segmentMismatch', ...
        '%s has %d segments, the other inputs have %d.', label, size(out, 2), Nsegs);
end
if size(out, 3) == 1 && Nmice > 1
    out = repmat(out, 1, 1, Nmice);   % one atlas volume, shared by every mouse
end
end

%--------------------------------------------------------------------------
function name = divisionOfTarget(hier, rows)
%DIVISIONOFTARGET The division a target sits in, for the third name column.
%   Mirrors REORGANIZEAREAS, whose third column is the coarsest level at every
%   aggregation depth.
name = '';
divnames = hier.division(rows);
divnames = divnames(~cellfun('isempty', divnames));
if ~isempty(divnames)
    name = divnames{1};
end
end
