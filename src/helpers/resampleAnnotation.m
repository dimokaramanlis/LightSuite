function [annnew, lut, regions] = resampleAnnotation(annvol, parcelinfo, targetLevel)
%RESAMPLEANNOTATION  Coarsen a CCF annotation volume up the parcellation hierarchy.
%
%   The Allen CCF annotation volume stores a fine "parcellation_index" in
%   every voxel. Each index belongs to a substructure, which in turn rolls
%   up through the ontology:
%
%       organ -> category -> division -> structure -> substructure
%      (   2  /     5     /    26     /    358     /      681      regions)
%
%   This function relabels ANNVOL so that every voxel carries a new, compact
%   index representing the region it belongs to at a COARSER level.
%
%   ANNNEW = RESAMPLEANNOTATION(ANNVOL, PARCELINFO) coarsens to 'structure'
%   (the level directly above substructure).
%
%   ANNNEW = RESAMPLEANNOTATION(ANNVOL, PARCELINFO, TARGETLEVEL) coarsens to
%   the requested level. TARGETLEVEL is one of (case-insensitive):
%       'organ' | 'category' | 'division' | 'structure' | 'substructure'
%
%   [ANNNEW, LUT, REGIONS] = RESAMPLEANNOTATION(...) also returns:
%       LUT     - the (maxIndex+1)-by-1 uint16 lookup vector, so the same
%                 mapping can be reapplied without recomputing:
%                     annnew = reshape(lut(double(othervol)+1), size(othervol));
%       REGIONS - a table (one row per new index) with the acronym, full
%                 name and RGB colour of each coarse region. Handy for
%                 building a colormap (see example below).
%
%   Inputs
%   ------
%   ANNVOL      Numeric array (2-D or 3-D) of parcellation_index values.
%   PARCELINFO  Table read from 'parcellation_to_parcellation_term_membership.csv',
%               e.g. parcelinfo = readtable(csvfile);
%
%   Conventions
%   -----------
%   * Voxels equal to 0, and any "unassigned" index, map to 0 in the output
%     (i.e. background / outside-brain stays background).
%   * Real regions are numbered 1..N in order of first appearance in the file.
%   * Indices in ANNVOL that are not present in PARCELINFO are mapped to 0.
%
%   Example
%   -------
%       parcelinfo = readtable('parcellation_to_parcellation_term_membership.csv');
%       av         = readNPY('annotation_10.npy');          % your fine volume
%       [avCoarse, ~, regions] = resampleAnnotation(av, parcelinfo, 'division');
%
%       % Build an RGB colormap indexed by the new labels (row 1 = background):
%       cmap = [0 0 0; regions.R regions.G regions.B] / 255;
%       rgb  = ind2rgb(avCoarse(:,:,300) + 1, cmap);        % colour one slice
%       imshow(rgb)

    % ------------------------------------------------------------------ %
    % 1. Parse / validate inputs
    % ------------------------------------------------------------------ %
    if nargin < 3 || isempty(targetLevel)
        targetLevel = 'structure';
    end
    targetLevel = char(string(targetLevel));

    required = {'parcellation_index','parcellation_term_set_name', ...
                'parcellation_term_label','parcellation_term_acronym', ...
                'parcellation_term_name'};
    missing = required(~ismember(required, parcelinfo.Properties.VariableNames));
    if ~isempty(missing)
        error('resampleAnnotation:missingColumns', ...
              'parcelinfo is missing required column(s): %s', strjoin(missing, ', '));
    end

    % ------------------------------------------------------------------ %
    % 2. Keep only the rows describing the requested hierarchy level
    % ------------------------------------------------------------------ %
    setNames = string(parcelinfo.parcellation_term_set_name);
    keep     = strcmpi(setNames, targetLevel);
    if ~any(keep)
        opts = unique(setNames, 'stable');
        error('resampleAnnotation:badLevel', ...
              'targetLevel "%s" not found. Available levels: %s', ...
              targetLevel, strjoin(cellstr(opts), ', '));
    end
    P = parcelinfo(keep, :);

    % ------------------------------------------------------------------ %
    % 3. Assign a new compact id to each distinct coarse region
    %    (drop "unassigned" -> it becomes background = 0)
    % ------------------------------------------------------------------ %
    isBG  = strcmpi(string(P.parcellation_term_acronym), 'unassigned');
    realP = P(~isBG, :);

    [uLabels, firstRow, grp] = unique(string(realP.parcellation_term_label), 'stable');
    nRegions = numel(uLabels);          %#ok<NASGU>  (documented for clarity)
    coarseId = grp;                     % 1..nRegions, one per real fine index

    % ------------------------------------------------------------------ %
    % 4. Build the lookup table indexed by (fineIndex + 1)
    %    (index 0 is never written, so lut(1) stays 0 = background)
    % ------------------------------------------------------------------ %
    fineIdx = double(realP.parcellation_index);
    maxIdx  = max(double(P.parcellation_index));        % covers the full index range
    lut     = zeros(maxIdx + 1, 1, 'uint16');
    lut(fineIdx + 1) = uint16(coarseId);

    % ------------------------------------------------------------------ %
    % 5. Apply the LUT to the whole volume in one vectorised shot
    % ------------------------------------------------------------------ %
    av        = double(annvol);
    outOfRange = ~isfinite(av) | av < 0 | av > maxIdx;  % stray / unknown indices
    av(outOfRange) = 0;
    annnew = reshape(lut(av + 1), size(annvol));

    % ------------------------------------------------------------------ %
    % 6. Region lookup table (new index -> acronym / name / colour)
    % ------------------------------------------------------------------ %
    if nargout > 2
        Rrep = realP(firstRow, :);      % one representative row per coarse region
        regions = table( (1:numel(uLabels)).', ...
                         string(Rrep.parcellation_term_acronym), ...
                         string(Rrep.parcellation_term_name), ...
                         double(Rrep.red), double(Rrep.green), double(Rrep.blue), ...
                         'VariableNames', {'index','acronym','name','R','G','B'});
    end
end