function [f, h] = plotSpinalCordGrid(cordGridData, ax, cord_figure, black_bg, show_labels)
%PLOTSPINALCORDGRID Draw the spinal cord wire mesh, the cord PLOTBRAINGRID.
%
%   PLOTSPINALCORDGRID() draws the outline of the spinal cord atlas in a new
%   figure: a ring per spinal segment with its gray matter inside, lines along
%   the cord tying the rings into a tube, and the name of every segment
%   underneath. The mesh is traced from the atlas by MAKESPINALCORDGRID.
%
%   PLOTSPINALCORDGRID(GRIDDATA, AX, FIG, BLACK_BG, SHOW_LABELS) uses a grid
%   made earlier by MAKESPINALCORDGRID, draws into the axes AX (or into a new
%   axes in the figure FIG), on a black background when BLACK_BG is true, and
%   without segment names when SHOW_LABELS is false. Any of them can be [].
%
%   The cord lies along the horizontal axis, rostral (C1) at the origin and
%   dorsal up, in micrometres. A point [x y z] in atlas voxels, as
%   TRANSFORMCORDPOINTSTOATLAS saves it, goes on the mesh with
%
%       res = GRIDDATA.atlasres;           % [10 10 20] um, per [row col plane]
%       plot3(z*res(3), x*res(2), y*res(1))
%
%   [F, H] = PLOTSPINALCORDGRID(...) returns the figure and a struct with the
%   handles of the outline, the gray matter and the labels.
%
%   See also MAKESPINALCORDGRID, PLOTBRAINGRID, VISUALIZECELLDETECTIONS.

if nargin < 1 || isempty(cordGridData)
    cordGridData = makeSpinalCordGrid();
end

if nargin < 2 || isempty(ax)
    if nargin < 3 || isempty(cord_figure)
        cord_figure = figure('Name', 'Spinal Cord View');
    end
    ax = axes('Parent', cord_figure);
end

black_bg    = nargin >= 4 && ~isempty(black_bg) && black_bg;
show_labels = nargin < 5 || isempty(show_labels) || show_labels;

if black_bg
    linecol = [.7 .7 .7];
    textcol = [.7 .7 .7];
    set(get(ax, 'Parent'), 'Color', 'k')
else
    linecol = [0 0 0];
    textcol = [.35 .35 .35];
end
%--------------------------------------------------------------------------
washeld = ishold(ax);
hold(ax, 'on');

% [x y z] of the atlas -> [z x y] on screen, so the cord lies along x
ol = cordGridData.outline;
gm = cordGridData.greymatter;
h  = struct();
h.outline    = plot3(ax, ol(:, 3), ol(:, 1), ol(:, 2), 'Color', [linecol 0.3]);
h.greymatter = plot3(ax, gm(:, 3), gm(:, 1), gm(:, 2), 'Color', [linecol 0.15]);

h.labels = gobjects(0);
if show_labels
    seg   = cordGridData.segments;
    Nseg  = height(seg);
    xmid  = (min(ol(:, 1)) + max(ol(:, 1))) / 2;
    ylab  = max(ol(:, 2)) + 0.15*(max(ol(:, 2)) - min(ol(:, 2)));
    h.labels = text(ax, seg.zring_um, repmat(xmid, Nseg, 1), repmat(ylab, Nseg, 1), ...
        cellstr(seg.name), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 7, 'Color', textcol);
end

if ~washeld
    hold(ax, 'off');
end
%--------------------------------------------------------------------------
% no 'vis3d' here, unlike the brain: it sizes the plot for any rotation, which
% for something 15 times longer than it is wide leaves a sliver of cord
set(ax, 'ZDir', 'reverse')
axis(ax, 'equal');
axis(ax, 'tight');
axis(ax, 'off');
view(ax, -10, 25);
f = get(ax, 'Parent');

end
