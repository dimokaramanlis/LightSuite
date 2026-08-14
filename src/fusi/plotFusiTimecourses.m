function cf = plotFusiTimecourses(lags, traces, opts)
%PLOTFUSITIMECOURSES Peristimulus block-triggered fUS responses (Fig. 1B style).
%   cf = plotFusiTimecourses(lags, traces) plots the mean +/- SEM
%   peristimulus response for the combined, object and scrambled stimulus
%   groups on one axis, with the 12 s stimulation window shaded.
%
%   Input:
%     lags    [nlags x 1] time from block onset (s).
%     traces  struct with fields combined/object/scrambled, each a struct with
%             .mean [nlags x 1] and optional .sem [nlags x 1]; OR a numeric
%             [nlags x 3] matrix of means in that column order.
%     opts    (optional): .stimdur (default 12), .title, .ylabel
%             (default '\DeltaI/I'), .savepng (''), .visible ('on').
%
%   See also FUSIPERISTIMULUSTIMECOURSES.

if nargin < 3, opts = struct(); end
opts = setdefault(opts, 'stimdur', 12);
opts = setdefault(opts, 'title', '');
opts = setdefault(opts, 'ylabel', '\DeltaI/I');
opts = setdefault(opts, 'savepng', '');
opts = setdefault(opts, 'visible', 'on');

lags   = lags(:);
groups = {'combined','object','scrambled'};
colors = [0 0 0; 0.85 0.33 0.10; 0.00 0.45 0.74];   % black / orange / blue

% normalize input to mean/sem per group
M = struct();
if isnumeric(traces)
    for g = 1:3, M.(groups{g}) = struct('mean', traces(:,g), 'sem', []); end
else
    M = traces;
end

% y-range from the data (mean +/- sem), padded, computed up front so the
% stimulation-window patch can be drawn once at the bottom of the stack
lo = inf; hi = -inf;
for g = 1:3
    if ~isfield(M, groups{g}), continue; end
    m = M.(groups{g}).mean(:); s = zeros(size(m));
    if isfield(M.(groups{g}), 'sem') && ~isempty(M.(groups{g}).sem), s = M.(groups{g}).sem(:); end
    lo = min(lo, min(m - s)); hi = max(hi, max(m + s));
end
if ~isfinite(lo), lo = -1; hi = 1; end
pad = 0.05 * (hi - lo + eps); yl = [lo - pad, hi + pad];

cf = figure('Visible', opts.visible, 'Color', 'w', 'Position', [60 60 720 460]);
ax = axes(cf); hold(ax, 'on');
% (1) stimulation window, drawn first so everything else sits on top
patch(ax, [0 opts.stimdur opts.stimdur 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
yline(ax, 0, 'k-', 'HandleVisibility', 'off');
% (2) SEM shading, then (3) mean lines
h = gobjects(1,3);
for g = 1:3
    if ~isfield(M, groups{g}), continue; end
    m = M.(groups{g}).mean(:); c = colors(g,:);
    if isfield(M.(groups{g}), 'sem') && ~isempty(M.(groups{g}).sem)
        s = M.(groups{g}).sem(:);
        patch(ax, [lags; flipud(lags)], [m-s; flipud(m+s)], c, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    h(g) = plot(ax, lags, m, 'Color', c, 'LineWidth', 1.8);
end
xlim(ax, [lags(1) lags(end)]); ylim(ax, yl);
xlabel(ax, 'time from block onset (s)'); ylabel(ax, opts.ylabel);
legend(ax, h(isgraphics(h)), groups(isgraphics(h)), 'Location', 'northeast', 'Box', 'off');
if ~isempty(opts.title), title(ax, opts.title, 'Interpreter', 'none'); end
box(ax, 'off');

if ~isempty(opts.savepng), print(cf, opts.savepng, '-dpng', '-r120'); end
end

% -------------------------------------------------------------------------
function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
