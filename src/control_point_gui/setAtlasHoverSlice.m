function setAtlasHoverSlice(htxt, avslice)
%SETATLASHOVERSLICE Give the hover label the annotation of the displayed slice.
%
%   SETATLASHOVERSLICE(HTXT, AVSLICE) stores AVSLICE (annotation values, same
%   size and orientation as the displayed atlas image) for the label made by
%   createAtlasHoverLabel, moves the label to the image's bottom right corner
%   and refreshes it, so the name is right after scrolling too.

if isempty(htxt) || ~isgraphics(htxt)
    return;
end
ud         = htxt.UserData;
ud.avslice = avslice;
ud.lastid  = nan;
htxt.UserData = ud;

ax = htxt.Parent;
% bottom right pixel edge of the image; text grows up and left from there
xr = size(avslice, 2) + 0.5;
if strcmp(ax.YDir, 'reverse')
    yb = size(avslice, 1) + 0.5;
else
    yb = 0.5;
end
htxt.Position = [xr, yb, 0];

refreshAtlasHoverLabel(htxt);
end
