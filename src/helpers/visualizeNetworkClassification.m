function fig = visualizeNetworkClassification(net, Xval, labels, pathsave)
%VISUALIZENETWORKCLASSIFICATION Show example cells/non-cells classified by a network.
%
%   visualizeNetworkClassification(NET, XVAL) runs the trained classification
%   network NET on the validation images XVAL and plots up to 100 examples
%   classified as "cells" and up to 100 classified as "non-cells". Each
%   detection is drawn as a grayscale composite of its views (side by side),
%   tiled into a montage in the style of the cell-classification figure panels.
%   The two montages are laid out with the panel library so the whitespace
%   gaps can be tuned via the margin parameters near the top of this file.
%
%   visualizeNetworkClassification(NET, XVAL, LABELS) additionally compares the
%   predictions against the ground-truth LABELS and outlines every mismatched
%   detection with a red box (false positives among the cells, false negatives
%   among the non-cells).
%
%   visualizeNetworkClassification(NET, XVAL, LABELS, PATHSAVE) also saves the
%   figure as a PNG inside the folder PATHSAVE. The file name carries the exact
%   date and time of creation, e.g. NetworkClassification_20260701_142530.png.
%   Pass LABELS as [] to skip mismatch marking but still save.
%
%   Inputs:
%     NET      - trained SeriesNetwork/DAGNetwork (labels 1 = cell, 0 = noise)
%     XVAL     - [H x W x C x Nval] image stack as produced by prepareImagesForCNN
%                (values already scaled to [0 1])
%     LABELS   - (optional) [Nval x 1] ground-truth labels (numeric 1/0, logical,
%                or categorical) used to mark mismatches
%     PATHSAVE - (optional) folder in which to save the PNG
%
%   Output:
%     FIG      - handle to the created figure
%
%   See also PREPAREIMAGESFORCNN, TRAINNETWORK, CLASSIFY, PANEL.

% -------------------------------------------------------------------------
% Parameters  (tweak whitespace / text / brightness here)
% -------------------------------------------------------------------------
maxPerClass = 100;          % max examples shown per class
nCols       = 10;           % montage columns
txtsize     = 7;            % font size for titles
whitePctile = 99;           % display white-point percentile of cell pixels
blackPctile  = 1;
                            %   (lower -> brighter cells, higher -> dimmer)
fw          = 32;           % figure width  (cm)
fh          = 22;           % figure height (cm)
outerMargin = [1 1 1 5];    % figure outer margin [L B R T] (mm)
panelGap    = 4;            % vertical gap between the two montages (mm)
dosave      = nargin > 3 && ~isempty(pathsave);
haveLabels  = nargin > 2 && ~isempty(labels);

% -------------------------------------------------------------------------
% Classify validation images (labels: 1 = cell, 0 = noise)
% -------------------------------------------------------------------------
YPred   = classify(net, Xval);
isCell  = (YPred == "1");
isNoise = (YPred == "0");
if ~any(isCell) && ~any(isNoise)
    % Fallback for unexpected class names: assume [noise cell] ordering
    cats    = categories(removecats(YPred));
    isNoise = (YPred == cats{1});
    if numel(cats) > 1
        isCell = (YPred == cats{end});
    end
end
nCellTot  = nnz(isCell);
nNoiseTot = nnz(isNoise);
nTot      = numel(YPred);

% Ground-truth "is a cell" per detection, for mismatch marking
if haveLabels
    trueIsCell = labelIsCell(labels);
    if numel(trueIsCell) ~= nTot
        warning('visualizeNetworkClassification:labelSize', ...
            'LABELS has %d entries but XVAL has %d images; skipping mismatch marking.', ...
            numel(trueIsCell), nTot);
        haveLabels = false;
    end
end

% -------------------------------------------------------------------------
% Subsample up to maxPerClass of each, reproducibly, without disturbing the
% caller's global RNG state
% -------------------------------------------------------------------------
oldRng = rng; rng(8);
idxCell  = pickSubset(find(isCell),  maxPerClass);
idxNoise = pickSubset(find(isNoise), maxPerClass);
rng(oldRng);

% Mismatch flags for the shown detections (in tile order): a shown "cell" is a
% mismatch when its label is noise (false positive); a shown "non-cell" is a
% mismatch when its label is cell (false negative).
if haveLabels
    mmCell  = ~trueIsCell(idxCell);
    mmNoise =  trueIsCell(idxNoise);
else
    mmCell  = [];
    mmNoise = [];
end

% -------------------------------------------------------------------------
% Build montages and pick a shared display range from the cell pixels so the
% cells are bright (non-cells are shown on the same scale, hence dimmer)
% -------------------------------------------------------------------------
[imgCell,  H, Wc] = buildMontage(Xval, idxCell,  nCols);
 imgNoise         = buildMontage(Xval, idxNoise, nCols);

vals = imgCell(imgCell > 0);
if isempty(vals), vals = imgNoise(imgNoise > 0); end
if isempty(vals)
    climUse = [0 0.1];
else
    climUse = prctile(vals, [blackPctile whitePctile]);
end

% -------------------------------------------------------------------------
% Plot with panel
% -------------------------------------------------------------------------
fig = figure('Color', 'w', 'Name', 'Network classification', ...
    'NumberTitle', 'off', 'Units', 'centimeters', 'Position', [2 2 fw fh]);
p = panel();
p.pack('v', 2);
p.fontsize  = txtsize;
p.margin    = outerMargin;
p.de.margin = panelGap;

p(1).select();
drawMontage(imgCell, H, Wc, climUse, mmCell);
title(classTitle('Cells', numel(idxCell), nCellTot, nTot, mmCell));

p(2).select();
drawMontage(imgNoise, H, Wc, climUse, mmNoise);
title(classTitle('Non-cells', numel(idxNoise), nNoiseTot, nTot, mmNoise));

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
if dosave
    if ~exist(pathsave, 'dir'); mkdir(pathsave); end
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    fout  = fullfile(pathsave, sprintf('NetworkClassification_%s.png', stamp));
    p.export(fout, sprintf('-w%d', fw * 10), sprintf('-h%d', fh * 10), '-r300');
    fprintf('Saved classification montage to %s\n', fout);
end

if nargout == 0
    clear fig
end
end

% =========================================================================
function s = classTitle(name, nShown, nTot, nAll, mm)
% Build a panel title, appending the mismatch count when labels are available.
s = sprintf('%s: %d shown of %d (%.1f%% of detections)', ...
    name, nShown, nTot, 100 * nTot / nAll);
if ~isempty(mm)
    s = sprintf('%s — %d mismatches (red)', s, nnz(mm));
end
end

% =========================================================================
function tf = labelIsCell(labels)
% Return a logical column 'is a cell (label 1)' for numeric / logical /
% categorical / string labels.
if iscategorical(labels)
    if ismember("1", string(categories(labels)))
        tf = (labels == "1");
    else
        cats = categories(removecats(labels));   % assume last category = cell
        tf   = (labels == cats{end});
    end
elseif islogical(labels)
    tf = labels;
else
    tf = (labels == 1);
end
tf = tf(:);
end

% =========================================================================
function sub = pickSubset(idx, nmax)
% Randomly keep at most NMAX entries of IDX.
if numel(idx) > nmax
    idx = idx(randperm(numel(idx), nmax));
end
sub = idx;
end

% =========================================================================
function [bigimg, H, Wc] = buildMontage(X, idx, nCols)
% Tile the detections in IDX into one image. Each tile is a detection's views
% concatenated side by side and normalized by its L2 norm (matching the
% manuscript figure panels). Empty tiles are left as NaN (drawn transparent).

[H, W, C, ~] = size(X);
Wc = W * C;                 % composite width (all views side by side)
n  = numel(idx);
if n == 0
    bigimg = [];
    return;
end

nCols  = min(nCols, n);
nRows  = ceil(n / nCols);
bigimg = nan(nRows * H, nCols * Wc, 'single');

for k = 1:n
    r    = floor((k - 1) / nCols);
    c    = mod(k - 1, nCols);
    comp = reshape(single(X(:, :, :, idx(k))), H, Wc);  % views side by side
    nrm  = sqrt(sum(comp(:).^2));
    if nrm > 0, comp = comp ./ nrm; end
    bigimg(r * H + (1:H), c * Wc + (1:Wc)) = comp;
end
end

% =========================================================================
function drawMontage(bigimg, H, Wc, climUse, mismatch)
% Draw a montage built by buildMontage on the current (panel) axes, outlining
% the tiles flagged in MISMATCH (logical, in tile order) with a red box.
if nargin < 5, mismatch = []; end
cla;
if isempty(bigimg)
    axis off;
    text(0.5, 0.5, 'No examples', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
    ax = gca; ax.Title.Visible = 'on';
    return;
end

nRows = size(bigimg, 1) / H;
nCols = size(bigimg, 2) / Wc;

him = imagesc(bigimg, climUse);
set(him, 'AlphaData', ~isnan(bigimg));   % empty tiles stay white
axis image off; hold on;
colormap(gca, gray);

% yellow grid between detections
xline(0.5 + (0:nCols) * Wc, 'Color', 'y');
yline(0.5 + (0:nRows) * H,  'Color', 'y');

% red box around each mismatched detection
for k = find(mismatch(:)')
    r = floor((k - 1) / nCols);
    c = mod(k - 1, nCols);
    rectangle('Position', [c * Wc + 1.5, r * H + 1.5, Wc - 2, H - 2], ...
        'EdgeColor', 'r', 'LineWidth', 1.5);
end
hold off;

ax = gca; ax.Title.Visible = 'on';       % keep title after 'axis off'
end
