function tbl = cordAreaTable(hier, segmentnames, values)
%CORDAREATABLE Long-format table of cord statistics, one row per region+segment.
%
%   TBL = CORDAREATABLE(HIER, SEGMENTNAMES, VALUES) builds the table every cord
%   CSV is written from. HIER comes from CORDAREAHIERARCHY and names each region
%   at all three levels; SEGMENTNAMES is the segment axis, GROUPING.SEGNAMES
%   from CORDATLASGROUPING; VALUES is a struct whose every field is an
%   Nareas x Nsegments array, its field name becoming a column.
%
%   The result carries the same identifying columns as the brain CSVs written by
%   GENERATEREGISTEREDBRAINVOLUMES - name, structure, division and
%   parcellation_index - so a region can be found by name rather than by looking
%   its id up in the atlas. Where the brain splits its measurements into a left
%   and a right column, the cord has a whole rostrocaudal axis to carry, so the
%   table is long instead of wide: one row per region and segment, with a
%   'segment' column. That keeps the file shape independent of how many segments
%   the atlas defines.
%
%   See also CORDAREAHIERARCHY, CORDATLASGROUPING, GENERATEREGISTEREDCORDVOLUME,
%   TRANSFORMCORDPOINTSTOATLAS.

Nareas    = height(hier);
Nsegments = numel(segmentnames);

[iarea, iseg] = ndgrid(1:Nareas, 1:Nsegments);
iarea = iarea(:);
iseg  = iseg(:);

tbl = table(hier.name(iarea), hier.acronym(iarea), ...
    hier.structure(iarea), hier.division(iarea), ...
    hier.parcellation_index(iarea), string(segmentnames(iseg)), ...
    'VariableNames', {'name', 'acronym', 'structure', 'division', ...
                      'parcellation_index', 'segment'});

valnames = fieldnames(values);
for ii = 1:numel(valnames)
    block = values.(valnames{ii});
    assert(isequal(size(block, 1), Nareas) && isequal(size(block, 2), Nsegments), ...
        'cordAreaTable:sizeMismatch', ...
        '%s is %s but the table needs %d areas x %d segments.', ...
        valnames{ii}, mat2str(size(block)), Nareas, Nsegments);
    tbl.(valnames{ii}) = double(block(:));
end

end
