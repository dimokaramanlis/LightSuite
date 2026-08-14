function cf = plotFusiActivationMontage(statvol, opts)
%PLOTFUSIACTIVATIONMONTAGE Coronal montage of a fUS activation map.
%   cf = plotFusiActivationMontage(statvol) shows the 3-D statistic volume
%   statvol (correlation, beta or T-score) as a grid of coronal slices, in
%   the spirit of Fig. 1C. With an anatomy underlay it draws the thresholded
%   statistic in colour on top of the greyscale anatomy.
%
%   opts (optional):
%     .underlay   [same size as statvol] anatomy for the greyscale backdrop
%                 (e.g. the session template or the per-mouse volumeavg). []
%                 for a plain coloured montage.
%     .clim       colour limits (default symmetric [-m m], m = 99th pctile).
%     .thresh     hide |stat| <= thresh over the underlay (default 0).
%     .sliceaxis  axis whose slices tile the montage (default 3, coronal).
%     .sliceidx   indices along sliceaxis to show (default all; use e.g.
%                 round(linspace(lo,hi,36)) for a large atlas volume).
%     .gridcols   montage columns (default auto).
%     .title      figure title.
%     .savepng    path (no extension) to print a PNG; '' to skip (default).
%     .visible    'on'/'off' (default 'on').
%
%   Returns the figure handle cf.
%
%   See also FUSISESSIONACTIVATIONMAPS, COMPUTEFUSIACTIVATIONMAPS.

if nargin < 2, opts = struct(); end
opts = setdefault(opts, 'underlay', []);
opts = setdefault(opts, 'thresh', 0);
opts = setdefault(opts, 'sliceaxis', 3);
opts = setdefault(opts, 'sliceidx', []);
opts = setdefault(opts, 'gridcols', []);
opts = setdefault(opts, 'title', '');
opts = setdefault(opts, 'savepng', '');
opts = setdefault(opts, 'visible', 'on');

statvol = single(statvol);
if ~isfield(opts, 'clim') || isempty(opts.clim)
    m = quantile(abs(statvol(isfinite(statvol))), 0.99);
    if ~isfinite(m) || m == 0, m = 1; end
    opts.clim = [-m m];
end

% orient so the chosen slice axis is the 3rd (montage) dimension
perm = 1:3; perm(perm == opts.sliceaxis) = []; perm = [perm opts.sliceaxis];
statvol = permute(statvol, perm);
if ~isempty(opts.underlay), opts.underlay = permute(opts.underlay, perm); end
% keep only the requested slices (e.g. to thin out a large atlas volume)
if ~isempty(opts.sliceidx)
    statvol = statvol(:, :, opts.sliceidx);
    if ~isempty(opts.underlay), opts.underlay = opts.underlay(:, :, opts.sliceidx); end
end
nsl     = size(statvol, 3);

cmap = divergingmap(256);                        % blue - white - red

% ---- build per-slice RGB ----------------------------------------------
[h, w, ~] = size(statvol);
rgb = zeros(h, w, 3, nsl, 'single');
if ~isempty(opts.underlay)
    U = single(opts.underlay);                 % already permuted + sliced
    ql = quantile(U(isfinite(U)), [0.01 0.99]);
    U  = min(max((U - ql(1)) / max(range(ql), eps), 0), 1);   % grey [0,1]
    U(~isfinite(U)) = 0;
end
for k = 1:nsl
    s   = statvol(:, :, k);
    idx = colorindex(s, opts.clim, size(cmap,1));
    col = ind2rgb(idx, cmap);                    % [h w 3] coloured stat
    if isempty(opts.underlay)
        rgb(:, :, :, k) = col;
    else
        g = repmat(U(:, :, k), [1 1 3]);
        a = single(abs(s) > opts.thresh);        % overlay where suprathresh
        a(~isfinite(s)) = 0;
        rgb(:, :, :, k) = g .* (1 - a) + col .* a;
    end
end

% ---- montage -----------------------------------------------------------
cf = figure('Visible', opts.visible, 'Color', 'w', 'Position', [50 50 1500 950]);
if isempty(opts.gridcols)
    montage(rgb);
else
    montage(rgb, 'Size', [ceil(nsl/opts.gridcols) opts.gridcols]);
end
colormap(gca, cmap); caxis(opts.clim); cb = colorbar; cb.Color = [0 0 0];
if ~isempty(opts.title), title(opts.title, 'Interpreter', 'none'); end

if ~isempty(opts.savepng)
    print(cf, opts.savepng, '-dpng', '-r120');
end
end

% -------------------------------------------------------------------------
function idx = colorindex(s, clim, n)
s = min(max(s, clim(1)), clim(2));
idx = round((s - clim(1)) / (clim(2) - clim(1)) * (n - 1)) + 1;
idx(~isfinite(idx)) = round(n/2);                % NaN -> mid (neutral) colour
end

function cmap = divergingmap(n)
%DIVERGINGMAP Blue -> white -> red, n rows.
h  = floor(n/2);
up = linspace(0, 1, h).';
lo = linspace(0, 1, n - h).';
cmap = [ [zeros(h,1)+0.1+0.9*up, zeros(h,1)+0.2+0.8*up, ones(h,1)]; ...   % blue->white
         [ones(n-h,1), 1-0.8*lo, 1-0.9*lo] ];                            % white->red
cmap = min(max(cmap, 0), 1);
end

function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
