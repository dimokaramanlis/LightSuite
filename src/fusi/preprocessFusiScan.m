function [X, info] = preprocessFusiScan(I, tframes, opts)
%PREPROCESSFUSISCAN Baseline, drift-removal and smoothing of a fUS recording.
%   [X, info] = preprocessFusiScan(I, tframes) converts a raw power-Doppler
%   recording into the voxels-by-time matrix used for activation mapping.
%   Steps (paper "fUS preprocessing"/"fUS activation", with one deliberate
%   change to the high-pass, see below):
%     1) relative change: dI/I = I ./ mean_t(I) - 1, per voxel, so voxels of
%        very different baseline Doppler amplitude become comparable.
%     2) drift removal: zero-phase Butterworth high-pass (filtfilt).
%     3) temporal smoothing: moving average over opts.smoothframes frames.
%
%   Input:
%     I        [nz nxy nt] or [nvox nt] raw power-Doppler (FUS.mat "I").
%     tframes  [nt x 1] frame times (FUS.mat "time"), used only for the rate.
%     opts     (optional) struct:
%        .hpcutoff    high-pass cutoff in Hz (default 0.01). NOTE: the paper
%                     states 0.056 Hz, but the block fundamental is 1/24 s =
%                     0.0417 Hz, which sits BELOW 0.056 Hz; a zero-phase
%                     filter at 0.056 Hz removes the response itself. 0.01 Hz
%                     removes slow drift while preserving the block band. Set
%                     to 0/[] to skip high-pass (linear detrend only).
%        .hporder     Butterworth order (default 2).
%        .smoothframes moving-average window in frames (default 4).
%        .baseline    'mean' (default) or 'median' temporal baseline.
%        .chunk       voxels filtered per batch (default 20000) to bound RAM.
%
%   Output:
%     X    [nt x nvox] single, preprocessed (time down columns; each column a
%          voxel), ready for corr / GLM. Layout matches computeFusiActivationMaps.
%     info struct with .fs, .b, .a (filter, [] if none), .volshape, .nt,
%          .baseline vector, .opts.
%
%   The filter coefficients info.b/info.a are returned so the design
%   regressors can be filtered identically before fitting (proper practice).
%
%   See also COMPUTEFUSIACTIVATIONMAPS, ALIGNFUSISTIMTOSCAN.

if nargin < 3, opts = struct(); end
opts = setdefault(opts, 'hpcutoff', 0.01);
opts = setdefault(opts, 'hporder', 2);
opts = setdefault(opts, 'smoothframes', 4);
opts = setdefault(opts, 'baseline', 'mean');
opts = setdefault(opts, 'chunk', 20000);

% ---- shape bookkeeping (time is always the LAST dimension) -------------
sz       = size(I);
nt       = sz(end);
volshape = sz(1:end-1);                       % [d1 d2 d3], [nz nxy], or [nvox]
if numel(volshape) < 2, volshape = []; end    % plain [nvox x nt] input
nvox     = prod(sz(1:end-1));
X        = single(reshape(I, nvox, nt)).';    % [nt x nvox]
clear I;

tframes = double(tframes(:));
fs      = 1 / median(diff(tframes));

% ---- (1) relative change dI/I -----------------------------------------
if strcmpi(opts.baseline, 'median')
    base = median(X, 1);
else
    base = mean(X, 1);
end
base(~isfinite(base) | base == 0) = NaN;     % avoid divide-by-zero blow-ups
X = X ./ base - 1;
X(~isfinite(X)) = 0;

% ---- (2) high-pass drift removal (chunked filtfilt) --------------------
b = []; a = [];
if ~isempty(opts.hpcutoff) && opts.hpcutoff > 0
    Wn = opts.hpcutoff / (fs/2);
    [b, a] = butter(opts.hporder, Wn, 'high');
    for c0 = 1:opts.chunk:nvox
        cc = c0:min(c0+opts.chunk-1, nvox);
        X(:, cc) = single(filtfilt(b, a, double(X(:, cc))));
    end
else
    X = single(detrend(double(X), 1));       % at least remove linear drift
end

% ---- (3) temporal smoothing -------------------------------------------
if opts.smoothframes > 1
    X = movmean(X, opts.smoothframes, 1);
end

info = struct('fs', fs, 'b', b, 'a', a, 'volshape', volshape, ...
    'nt', nt, 'baseline', base, 'opts', opts);
end

% -------------------------------------------------------------------------
function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
