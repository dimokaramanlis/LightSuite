function markers = readCellCounterXML(xmlpath, coordoffset)
%READCELLCOUNTERXML Read an ImageJ Cell Counter marker file into point arrays.
%   MARKERS = READCELLCOUNTERXML(XMLPATH) parses a CellCounter_Marker_File XML
%   (the output of the ImageJ / Fiji "Cell Counter" plugin) and returns one
%   entry per <Marker_Type> block found in it.
%
%   A Cell Counter file looks like this:
%
%       <CellCounter_Marker_File>
%         <Image_Properties> ... </Image_Properties>
%         <Marker_Data>
%           <Current_Type>1</Current_Type>
%           <Marker_Type>
%             <Type>1</Type>
%             <Marker><MarkerX>2565</MarkerX><MarkerY>3734</MarkerY><MarkerZ>62</MarkerZ></Marker>
%             ...
%           </Marker_Type>
%           <Marker_Type> <Type>2</Type> ... </Marker_Type>
%         </Marker_Data>
%       </CellCounter_Marker_File>
%
%   Marker types are the plugin's counter categories. They are NOT image
%   channels: nothing in the file records which channel a type was counted on.
%   transformPointsToAtlas therefore maps type -> channel explicitly and prints
%   the mapping it used; see its 'channel' option.
%
%   COORDINATE CONVENTION. The plugin stores MarkerX/MarkerY as 0-based canvas
%   pixels and MarkerZ as the 1-based ImageJ slice number. LightSuite's
%   cell_locations arrays are 1-based in all three axes, so COORDOFFSET is
%   added to the raw values; it defaults to [1 1 0], which converts the plugin's
%   convention to LightSuite's. Pass [0 0 0] to take the file's numbers as they
%   are (e.g. for markers written by another tool).
%
%   Input:
%     xmlpath      path to the .xml file.
%     coordoffset  (optional) 1x3 [dx dy dz] added to every marker.
%                  Default [1 1 0].
%
%   Output:
%     markers      struct array, one element per <Marker_Type>:
%        .type     the <Type> value of the block (numeric)
%        .points   [N x 3] single, [x y z] in the cell_locations convention
%        .npoints  N
%
%   Empty <Marker_Type> blocks (a category the user never clicked) are dropped.
%
%   The file is parsed with sscanf over the whole marker block, which reads the
%   ~1.3 million markers of a large file in a few seconds. A slower regular
%   expression is used as a fallback if a block does not match the plugin's
%   usual layout.
%
%   See also TRANSFORMPOINTSTOATLAS, READCELLLOCATIONSFILE.

if nargin < 2 || isempty(coordoffset)
    coordoffset = [1 1 0];
end
coordoffset = double(coordoffset(:)).';
if numel(coordoffset) ~= 3
    error('readCellCounterXML:badOffset', 'coordoffset must be a 1x3 vector.');
end

if ~isfile(xmlpath)
    error('readCellCounterXML:fileNotFound', 'XML file not found: %s', xmlpath);
end

txt = fileread(xmlpath);

blockstarts = strfind(txt, '<Marker_Type>');
blockstops  = strfind(txt, '</Marker_Type>');
if isempty(blockstarts)
    error('readCellCounterXML:noMarkerTypes', ...
        ['No <Marker_Type> blocks found in %s. Is this an ImageJ Cell Counter ' ...
         'marker file?'], xmlpath);
end
if numel(blockstops) ~= numel(blockstarts)
    error('readCellCounterXML:malformed', ...
        'Found %d <Marker_Type> but %d </Marker_Type> in %s.', ...
        numel(blockstarts), numel(blockstops), xmlpath);
end

markers = struct('type', {}, 'points', {}, 'npoints', {});
for k = 1:numel(blockstarts)
    blk = txt(blockstarts(k):blockstops(k) - 1);

    % the block's own <Type> element, and where the markers start after it
    typetok = regexp(blk, '<Type>\s*(-?\d+)\s*</Type>', 'tokens', 'once');
    typeend = strfind(blk, '</Type>');
    if isempty(typetok) || isempty(typeend)
        warning('readCellCounterXML:noType', ...
            'Marker_Type block %d has no <Type>; numbering it %d.', k, k);
        thistype = k;
        body     = blk;
    else
        thistype = str2double(typetok{1});
        body     = blk(typeend(1) + numel('</Type>'):end);
    end

    nexpected = numel(strfind(body, '<Marker>'));
    if nexpected == 0
        continue;   % a counter category the user never used
    end

    pts = parseMarkerBody(body, nexpected, xmlpath, thistype);

    markers(end+1) = struct(...
        'type',    thistype, ...
        'points',  single(pts + coordoffset), ...
        'npoints', size(pts, 1)); %#ok<AGROW>
end

if isempty(markers)
    error('readCellCounterXML:noMarkers', ...
        'No markers found in any <Marker_Type> block of %s.', xmlpath);
end

end

%==========================================================================
function pts = parseMarkerBody(body, nexpected, xmlpath, thistype)
%PARSEMARKERBODY [N x 3] marker coordinates from one <Marker_Type> body.
%   Fast path: sscanf against the plugin's fixed element order. Whitespace in
%   the format matches any run of whitespace, so indentation and line breaks
%   do not matter, but the element order does.

fmt = [' <Marker> <MarkerX>%f</MarkerX> <MarkerY>%f</MarkerY>' ...
       ' <MarkerZ>%f</MarkerZ> </Marker>'];
v   = sscanf(body, fmt, [3 Inf]);

if size(v, 2) == nexpected
    pts = double(v.');
    return;
end

% Fallback: the block is laid out differently (extra elements, reordered
% axes, attributes). Pull each axis out by name instead - much slower, but it
% does not care about ordering.
warning('readCellCounterXML:slowParse', ...
    ['Marker_Type %d of %s does not follow the usual X/Y/Z element order ' ...
     '(fast parse read %d of %d markers); falling back to a regular ' ...
     'expression.'], thistype, xmlpath, size(v, 2), nexpected);

xs = regexp(body, '<MarkerX>\s*(-?[\d.eE+]+)\s*</MarkerX>', 'tokens');
ys = regexp(body, '<MarkerY>\s*(-?[\d.eE+]+)\s*</MarkerY>', 'tokens');
zs = regexp(body, '<MarkerZ>\s*(-?[\d.eE+]+)\s*</MarkerZ>', 'tokens');

if ~isequal(numel(xs), numel(ys), numel(zs))
    error('readCellCounterXML:axisMismatch', ...
        ['Marker_Type %d of %s has %d X, %d Y and %d Z entries; the file ' ...
         'looks truncated or malformed.'], ...
        thistype, xmlpath, numel(xs), numel(ys), numel(zs));
end

pts = [str2double([xs{:}]).', str2double([ys{:}]).', str2double([zs{:}]).'];

end
