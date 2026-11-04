function refreshAtlasHoverLabel(htxt)
%REFRESHATLASHOVERLABEL Update the hover label from the current mouse position.
%   Called on mouse motion (see createAtlasHoverLabel) and when the atlas
%   slice changes (see setAtlasHoverSlice).

if isempty(htxt) || ~isgraphics(htxt)
    return;
end
ud = htxt.UserData;
if isempty(ud.avslice)
    return;
end

cp  = htxt.Parent.CurrentPoint;
col = round(cp(1, 1));
row = round(cp(1, 2));
[nr, nc] = size(ud.avslice);

id = 0;
if row >= 1 && row <= nr && col >= 1 && col <= nc
    id = double(ud.avslice(row, col));
end
if isequal(id, ud.lastid)
    return;   % same region, nothing to redraw
end
ud.lastid     = id;
htxt.UserData = ud;

name = "";
if id ~= 0
    name = ud.lut.names(ud.lut.ids == id);
    if isempty(name)
        name = sprintf("id %d", id);
    end
end
if strlength(name) == 0
    set(htxt, 'String', '', 'Visible', 'off');
else
    set(htxt, 'String', char(name(1)), 'Visible', 'on');
end
end
