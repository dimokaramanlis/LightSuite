function [hier, levels] = cordAreaHierarchy(parcelinfo, avinds)
%CORDAREAHIERARCHY Three-level region hierarchy of the spinal cord atlas.
%
%   HIER = CORDAREAHIERARCHY(PARCELINFO, AVINDS) places every region id in
%   AVINDS - the finest grain the annotation volume actually labels - inside the
%   same three-level hierarchy the Allen brain parcellation uses, so cord and
%   brain outputs can be read the same way:
%
%     substructure : the labelled region itself, e.g. lamina 4 of the dorsal horn
%     structure    : laminae I-X, each combining its per-horn parts, and the
%                    dorsal (df), lateral (lf) and ventral (vf) funiculi
%     division     : gray matter (GM) and white matter (WM)
%
%   PARCELINFO is the Atlas_Regions.csv table of the Fiederling et al. (2021)
%   atlas (see LOADSPINALCORDATLASTABLES). AVINDS is the list of ids present in
%   the annotation, i.e. GROUPING.AVINDS from CORDATLASGROUPING.
%
%   The hierarchy is an *assignment*, not a set of overlapping queries: every id
%   gets at most one structure and one division, so summing over a level can
%   never count a region twice. Ids the atlas leaves outside the hierarchy - id
%   0 above all, which is everything outside the cord - keep empty structure and
%   division entries and are simply skipped when a level is aggregated.
%
%   HIER is a table with one row per entry of AVINDS, in the order given:
%       parcellation_index  atlas id of the substructure
%       acronym, name       its acronym and name
%       structure_index, structure_acronym, structure   its structure
%       division_index, division_acronym, division      its division
%
%   [HIER, LEVELS] = CORDAREAHIERARCHY(...) also returns the structure and
%   division targets themselves, as LEVELS.structure and LEVELS.division tables
%   with columns index/acronym/name. They come in anatomical order - laminae
%   I-X, then df, lf, vf; GM, then WM - which is the order the aggregated
%   outputs are written in, rather than in order of atlas id.
%
%   Nothing here depends on the spinal segment: a segment is a range of planes,
%   not a label in the annotation, so a region and a segment are independent
%   coordinates of the same measurement (see CORDATLASGROUPING). The hierarchy
%   therefore applies unchanged whatever the segment count of the atlas.
%
%   See also REORGANIZESPINALCORDAREAS, CORDATLASGROUPING,
%   LOADSPINALCORDATLASTABLES, REORGANIZEAREAS.

%--------------------------------------------------------------------------
% 1. the atlas table, normalised
%--------------------------------------------------------------------------
if ~ismember('id', parcelinfo.Properties.VariableNames)
    error('cordAreaHierarchy:noIdColumn', ...
        'parcelinfo must have an ''id'' column; it does not look like Atlas_Regions.csv.');
end

ids      = double(parcelinfo.id(:));
acronyms = textColumn(parcelinfo, 'acronym', numel(ids));
names    = textColumn(parcelinfo, 'name',    numel(ids));

% id -> row, the same trick cordAtlasGrouping uses for the annotation
rowlut          = zeros(max([ids; 0]) + 1, 1);
rowlut(ids + 1) = 1:numel(ids);

kids = childrenOfEveryRow(parcelinfo, ids, rowlut);

%--------------------------------------------------------------------------
% 2. the targets of each level
%--------------------------------------------------------------------------
% The cord hierarchy is curated, not derived: the direct children of GM are the
% two horns, while the level we want under it is the laminae, which the atlas
% expresses as the "combined" regions 201-210 listing their parts in
% children_IDs. White matter is the plain case, its children being the funiculi.
laminakeys   = arrayfun(@(x) {x}, 201:210, 'UniformOutput', false);
funiculikeys = {{'df'}, {'lf', 'lfc'}, {'vf'}};

levels.structure = resolveTargets([laminakeys, funiculikeys], ids, acronyms, names, rowlut);
levels.division  = resolveTargets({{'GM'}, {'WM'}},            ids, acronyms, names, rowlut);

%--------------------------------------------------------------------------
% 3. assign every labelled id to one structure and one division
%--------------------------------------------------------------------------
avinds = double(avinds(:));
Nav    = numel(avinds);
blanks = repmat({''}, Nav, 1);

hier = table(avinds, blanks, blanks, ...
    nan(Nav, 1), blanks, blanks, ...
    nan(Nav, 1), blanks, blanks, ...
    'VariableNames', {'parcellation_index', 'acronym', 'name', ...
        'structure_index', 'structure_acronym', 'structure', ...
        'division_index',  'division_acronym',  'division'});

for ii = 1:Nav
    row = rowOfId(avinds(ii), rowlut);
    if row > 0
        hier.acronym{ii} = acronyms{row};
        hier.name{ii}    = names{row};
    end
end

% id 0 is not a region, it is everything the annotation does not label
isbackground = avinds == 0;
hier.name(isbackground)    = {'outside cord'};
hier.acronym(isbackground) = {'bg'};

hier = assignLevel(hier, levels.structure, kids, rowlut, ...
    {'structure_index', 'structure_acronym', 'structure'});
hier = assignLevel(hier, levels.division,  kids, rowlut, ...
    {'division_index',  'division_acronym',  'division'});

%--------------------------------------------------------------------------
% 4. say what fell through, once, rather than per consumer
%--------------------------------------------------------------------------
% aggregating several channels calls this repeatedly with the same atlas, and
% the same warning ten times over is noise rather than information
persistent warnedabout
orphans = ~isbackground & isnan(hier.division_index);
orphanids = sprintf('%d,', hier.parcellation_index(orphans));
if any(orphans) && ~isequal(warnedabout, orphanids)
    warnedabout = orphanids;
    warning('cordAreaHierarchy:unassignedRegions', ...
        ['%d of %d labelled regions sit outside the gray/white matter ' ...
         'hierarchy and are dropped when aggregating (ids %s). This is ' ...
         'expected for the central canal and similar.'], ...
        nnz(orphans), Nav - nnz(isbackground), orphanids(1:end-1));
end

end

%==========================================================================
% Local helpers
%==========================================================================
function col = textColumn(tbl, varname, nrows)
%TEXTCOLUMN One column of a table as a cellstr, empty when the column is absent.
if ~ismember(varname, tbl.Properties.VariableNames)
    col = repmat({''}, nrows, 1);
    return
end
col = tbl.(varname);
if isstring(col) || ischar(col) || isnumeric(col)
    col = cellstr(string(col));
elseif ~iscell(col)
    col = repmat({''}, nrows, 1);
end
col = col(:);
col(~cellfun(@ischar, col)) = {''};
end

%--------------------------------------------------------------------------
function kids = childrenOfEveryRow(parcelinfo, ids, rowlut)
%CHILDRENOFEVERYROW Direct children of each row, from both linkage columns.
%   The atlas states parentage twice and not redundantly: parent_ID gives the
%   anatomical tree (horn -> lamina), while children_IDs is how the "combined"
%   regions name the parts they gather. Both are needed.

Nrow = numel(ids);
kids = cell(Nrow, 1);

if ismember('parent_ID', parcelinfo.Properties.VariableNames)
    parents = double(parcelinfo.parent_ID(:));
    for ii = 1:Nrow
        prow = rowOfId(parents(ii), rowlut);
        if prow > 0 && prow ~= ii
            kids{prow}(end + 1) = ids(ii);
        end
    end
end

if ismember('children_IDs', parcelinfo.Properties.VariableNames)
    childcol = parcelinfo.children_IDs;
    for ii = 1:Nrow
        listed = parseIdList(childcol, ii);
        listed(listed == ids(ii)) = [];   % a region is not its own child
        kids{ii} = [kids{ii}, listed];
    end
end

kids = cellfun(@unique, kids, 'UniformOutput', false);
end

%--------------------------------------------------------------------------
function out = parseIdList(col, irow)
%PARSEIDLIST The ids in one children_IDs cell, whatever readtable made of it.
if isnumeric(col)
    out = double(col(irow));
elseif iscell(col)
    out = idsOfValue(col{irow});
else
    out = idsOfValue(col(irow));
end
out = out(:)';
out = out(~isnan(out));
end

%--------------------------------------------------------------------------
function out = idsOfValue(val)
%IDSOFVALUE The numbers inside one children_IDs entry.
if isnumeric(val)
    out = double(val);
elseif ischar(val) || isstring(val)
    out = str2double(regexp(char(val), '\d+', 'match'));
else
    out = [];
end
end

%--------------------------------------------------------------------------
function row = rowOfId(id, rowlut)
%ROWOFID Row of an atlas id, 0 when the id is not in the table.
row = 0;
if isempty(id) || isnan(id) || id < 0 || id + 1 > numel(rowlut)
    return
end
row = rowlut(id + 1);
end

%--------------------------------------------------------------------------
function targets = resolveTargets(keys, ids, acronyms, names, rowlut)
%RESOLVETARGETS Turn the curated key list into the rows the atlas actually has.
%   Each key is a cell of alternatives, tried in order, so an atlas that spells
%   the lateral funiculus 'lfc' still resolves. A key with no match is dropped
%   with a warning rather than becoming an empty output row.

index   = nan(numel(keys), 1);
acronym = repmat({''}, numel(keys), 1);
name    = repmat({''}, numel(keys), 1);
missing = strings(0, 1);

for ii = 1:numel(keys)
    row = 0;
    for jj = 1:numel(keys{ii})
        alt = keys{ii}{jj};
        if isnumeric(alt)
            row = rowOfId(alt, rowlut);
        else
            hit = find(strcmpi(acronyms, alt), 1);
            if ~isempty(hit); row = hit; end
        end
        if row > 0; break; end
    end

    if row == 0
        missing(end + 1) = string(keys{ii}{1}); %#ok<AGROW>
        continue
    end

    index(ii)   = ids(row);
    acronym{ii} = acronyms{row};
    name{ii}    = names{row};
end

persistent warnedmissing
if ~isempty(missing) && ~isequal(warnedmissing, missing)
    warnedmissing = missing;
    warning('cordAreaHierarchy:missingTargets', ...
        ['%d hierarchy targets are not in the atlas table and are left out ' ...
         '(%s). Check that Atlas_Regions.csv is the Fiederling one.'], ...
        numel(missing), strjoin(missing, ', '));
end

keep    = ~isnan(index);
targets = table(index(keep), acronym(keep), name(keep), ...
    'VariableNames', {'index', 'acronym', 'name'});
end

%--------------------------------------------------------------------------
function hier = assignLevel(hier, targets, kids, rowlut, outvars)
%ASSIGNLEVEL Attach each substructure to the first target that contains it.
%   Walking down from the targets - rather than up from each id - keeps the
%   assignment single-valued by construction: the first target claiming an id
%   keeps it, so nothing is aggregated twice.

for ii = 1:height(targets)
    members = [targets.index(ii); descendantsOf(targets.index(ii), kids, rowlut)];

    take = ismember(hier.parcellation_index, members) & isnan(hier.(outvars{1}));
    if ~any(take); continue; end

    hier.(outvars{1})(take) = targets.index(ii);
    hier.(outvars{2})(take) = targets.acronym(ii);
    hier.(outvars{3})(take) = targets.name(ii);
end
end

%--------------------------------------------------------------------------
function out = descendantsOf(rootid, kids, rowlut)
%DESCENDANTSOF Every id below ROOTID, breadth first.
%   Iterative and visited-guarded: the atlas table contains rows that list
%   themselves among their children, which a plain recursion never returns from.

seen  = rootid;
queue = rootid;
while ~isempty(queue)
    current  = queue(1);
    queue(1) = [];
    row = rowOfId(current, rowlut);
    if row == 0; continue; end
    fresh = setdiff(kids{row}, seen);
    seen  = [seen, fresh];   %#ok<AGROW>
    queue = [queue, fresh];  %#ok<AGROW>
end
out = setdiff(seen(:), rootid);
end
