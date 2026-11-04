function lut = atlasNameLookup(pinfo, level)
%ATLASNAMELOOKUP Annotation value -> region name, for any supported atlas.
%
%   LUT = ATLASNAMELOOKUP(PINFO) builds a lookup from the region table
%   returned by loadAtlasInfo. LUT has fields
%       ids    - annotation values (double column)
%       names  - region names (string column), aligned with ids
%   Look a value up with  LUT.names(LUT.ids == value).
%
%   PINFO can be
%     * the Allen parcellation_to_parcellation_term_membership.csv table;
%       LEVEL picks the term set (default 'substructure', the level stored
%       in annotation_10)
%     * the spinal cord Atlas_Regions.csv table (columns id, name)
%     * a table with columns index and name, e.g. the REGIONS output of
%       resampleAnnotation (coarsened annotation ids)

if nargin < 2 || isempty(level)
    level = 'substructure';
end
vars = pinfo.Properties.VariableNames;

if ismember('parcellation_index', vars)
    % Allen: one row per (index, term set)
    keep  = strcmpi(string(pinfo.parcellation_term_set_name), level);
    ids   = double(pinfo.parcellation_index(keep));
    names = string(pinfo.parcellation_term_name(keep));
elseif all(ismember({'id', 'name'}, vars))
    % spinal cord atlas
    ids   = double(pinfo.id);
    names = string(pinfo.name);
elseif all(ismember({'index', 'name'}, vars))
    ids   = double(pinfo.index);
    names = string(pinfo.name);
else
    error('atlasNameLookup:unknownTable', ...
        'Region table not recognised (expected Allen, spinal cord or index/name columns).');
end

[ids, iu] = unique(ids(:));
names     = names(iu);
names(ismissing(names)) = "";
lut = struct('ids', ids, 'names', names(:));
end
