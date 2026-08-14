function tc = fusiPeristimulusTimecourses(X, tframes, blocks, offset, opts)
%FUSIPERISTIMULUSTIMECOURSES Block-triggered (peristimulus) fUS responses.
%   tc = fusiPeristimulusTimecourses(X, tframes, blocks, offset) averages the
%   preprocessed fUS signal in a window around each stimulus block onset
%   ("wrapped around presentations"), separately for the three groups used in
%   Fig. 1B: combined (all blocks), object and scrambled. It returns both a
%   per-voxel peristimulus volume-timecourse (which can be warped to the atlas
%   like the maps) and, if an ROI is given, a mean +/- SEM trace across
%   presentations.
%
%   Each block onset is placed on the scan clock at onset + offset (the same
%   offset used for the regressors), and the signal is sampled on a common
%   lag grid by linear interpolation, so irregular frame timing is handled
%   exactly. By default each epoch is baselined to its pre-onset mean.
%
%   Input:
%     X        [nt x nvox] preprocessed data (preprocessFusiScan).
%     tframes  [nt x 1] frame times (FUS.mat "time").
%     blocks   struct from loadFusiStimBlocks.
%     offset   stimulus->scan offset (s), e.g. -blocks.triggerPsyTime.
%     opts     (optional):
%        .window   [tpre tpost] s around onset (default [-6 24]; block is 12 s).
%        .dt       lag spacing s (default median frame period).
%        .baseline subtract per-epoch pre-onset mean (default true).
%        .roimask  logical [nvox x 1] ROI for the summary trace (default:
%                  none -> tc.roi is empty).
%
%   Output struct tc:
%     .lags       [nlags x 1] time from block onset (s)
%     .perVox     struct .combined/.object/.scrambled: [nvox x nlags] mean
%                 peristimulus response per voxel (reshape to volume with
%                 blocks-> the map volshape as needed).
%     .roi        struct (empty if no roimask) .combined/.object/.scrambled,
%                 each with .mean/.sem [nlags x 1] and .n (#presentations).
%     .nblocks    [combined object scrambled] presentation counts.
%     .window, .offset
%
%   See also LOADFUSISTIMBLOCKS, PREPROCESSFUSISCAN, FUSISESSIONACTIVATIONMAPS.

if nargin < 5, opts = struct(); end
opts = setdefault(opts, 'window',  [-6 24]);
opts = setdefault(opts, 'dt',      median(diff(tframes)));
opts = setdefault(opts, 'baseline', true);
opts = setdefault(opts, 'roimask', []);

tframes = double(tframes(:));
lags    = (opts.window(1):opts.dt:opts.window(2)).';
nlags   = numel(lags);
nvox    = size(X, 2);
baseIdx = lags < 0;                          % pre-onset baseline samples
haveROI = ~isempty(opts.roimask);

onset = blocks.onset(:) + offset;            % block onset on the scan clock
isObj = blocks.isObject(:);
isScr = blocks.isScrambled(:);

% per-voxel accumulators (combined = object + scrambled)
oSum = zeros(nlags, nvox, 'single'); oCnt = zeros(nlags, 1);
sSum = zeros(nlags, nvox, 'single'); sCnt = zeros(nlags, 1);
% per-block ROI traces (kept for mean/SEM across presentations)
if haveROI
    roiO = nan(nlags, sum(isObj)); roiS = nan(nlags, sum(isScr));
    jo = 0; js = 0;
end

for ib = 1:numel(onset)
    q     = onset(ib) + lags;                            % query times (scan)
    valid = q >= tframes(1) & q <= tframes(end);         % in-recording lags
    E     = interp1(tframes, X, q, 'linear');            % [nlags x nvox]
    if opts.baseline
        bl = baseIdx & valid;
        if any(bl), E = E - mean(E(bl, :), 1, 'omitnan'); end
    end
    Ev = E; Ev(~valid, :) = 0;                           % zero the clipped lags
    if isObj(ib)
        oSum = oSum + Ev; oCnt = oCnt + valid;
        if haveROI, jo = jo + 1; roiO(:, jo) = roimean(E, opts.roimask); end
    else
        sSum = sSum + Ev; sCnt = sCnt + valid;
        if haveROI, js = js + 1; roiS(:, js) = roimean(E, opts.roimask); end
    end
end

% per-voxel means (combined merges object + scrambled)
tc.perVox.object    = oSum ./ max(oCnt, 1);
tc.perVox.scrambled = sSum ./ max(sCnt, 1);
tc.perVox.combined  = (oSum + sSum) ./ max(oCnt + sCnt, 1);
% orient as [nvox x nlags] to match the map layout
tc.perVox.object    = tc.perVox.object.';
tc.perVox.scrambled = tc.perVox.scrambled.';
tc.perVox.combined  = tc.perVox.combined.';

tc.roi = struct([]);
if haveROI
    tc.roi = struct('object',    tracestat(roiO), ...
                    'scrambled', tracestat(roiS), ...
                    'combined',  tracestat([roiO roiS]));
end
tc.lags    = lags;
tc.nblocks = [sum(isObj)+sum(isScr), sum(isObj), sum(isScr)];
tc.window  = opts.window;
tc.offset  = offset;
end

% -------------------------------------------------------------------------
function m = roimean(E, mask)
m = mean(E(:, mask), 2, 'omitnan');
end

function s = tracestat(R)
n = sum(~isnan(R), 2);
s = struct('mean', mean(R, 2, 'omitnan'), ...
           'sem',  std(R, 0, 2, 'omitnan') ./ sqrt(max(n, 1)), ...
           'n',    size(R, 2));
end

function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
