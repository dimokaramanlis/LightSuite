function htxt = createAtlasHoverLabel(ax, lut)
%CREATEATLASHOVERLABEL Show the atlas region under the mouse in an atlas panel.
%
%   HTXT = CREATEATLASHOVERLABEL(AX, LUT) adds a text label to the image axes
%   AX, at the bottom right corner of the image (right-aligned to its edge),
%   and installs a WindowButtonMotionFcn on the figure that updates it with
%   the region name under the mouse. LUT comes from atlasNameLookup.
%
%   The GUI must hand over the annotation slice that matches the displayed
%   atlas image whenever that image changes:
%       setAtlasHoverSlice(HTXT, AVSLICE)
%   AVSLICE has the size of the displayed image, in the same orientation.

fig  = ancestor(ax, 'figure');
htxt = text(ax, nan, nan, '', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'Color', 'w', 'BackgroundColor', 'k', 'Margin', 2, ...
    'FontSize', 10, 'Interpreter', 'none', ...
    'HitTest', 'off', 'PickableParts', 'none', 'Clipping', 'off');
htxt.UserData = struct('lut', lut, 'avslice', [], 'lastid', nan);

set(fig, 'WindowButtonMotionFcn', @(~, ~) refreshAtlasHoverLabel(htxt));
end
