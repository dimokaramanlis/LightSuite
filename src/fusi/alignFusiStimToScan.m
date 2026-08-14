function [offset, diag] = alignFusiStimToScan(X, tframes, blocks, opts)
%ALIGNFUSISTIMTOSCAN Estimate the stimulus->scan clock offset for a session.
%   [offset, diag] = alignFusiStimToScan(X, tframes, blocks) finds the time
%   offset (seconds) that places the PsychoPy stimulus blocks onto the fUS
%   scan clock. It is needed because the PsychoPy log and the fUS acquisition
%   do NOT share a clock origin: the true offset is a per-session constant
%   (e.g. -60 s for DS_WT61/220221) and getting it wrong yields noise maps.
%
%   Method: for a grid of candidate offsets, build the HRF-convolved COMBINED
%   block regressor, correlate it with every voxel of the preprocessed data
%   X, and take the maximum voxel correlation as the alignment score. The
%   combined regressor is 24 s-periodic (block cadence), so its score peaks
%   repeat every 24 s; the absolute phase is then disambiguated with the
%   OBJECT-only regressor, whose object/scrambled sequence is pseudorandom
%   and therefore non-periodic. Finally the chosen offset is refined on a
%   fine grid. offset is defined so scanTime = psyTime + offset.
%
%   Input:
%     X        [nt x nvox] preprocessed data from preprocessFusiScan.
%     tframes  [nt x 1] frame times (FUS.mat "time").
%     blocks   struct from loadFusiStimBlocks.
%     opts     (optional):
%        .b,.a         regressor high-pass filter (match the data; from
%                      preprocessFusiScan info). Default [] (no filtering).
%        .hrfparams    HRF params (default [1.5 10 0.5 1 20 0 16]).
%        .coarsestep   coarse scan step, s (default 1).
%        .finestep     refine step, s (default 0.25).
%        .finehalfwin  refine half-window, s (default 3).
%        .searchwin    [lo hi] offsets to scan, s. Default: the feasible
%                      range in which all blocks still fall inside the scan,
%                      widened by 5 s.
%        .prioroffset  offset prior (s). Best: the key-press trigger,
%                      -blocks.triggerPsyTime (~1 s). Coarser: the wall-clock
%                      (psyStartDatenum - time0)*86400 (~+/-30 s). When finite
%                      it centres the search (+/- priorhalfwin) and is
%                      cross-checked against the result. Default NaN (unused).
%        .priorhalfwin half-width around prioroffset to scan. Use a small
%                      value (~10 s) for the key-press prior, ~90 s for the
%                      coarse wall-clock prior. Keep it < 12 s of the true
%                      value to avoid the 24 s block aliasing; within a wider
%                      window the object regressor resolves the alias.
%        .verbose      print progress (default true).
%
%   Output:
%     offset   scalar best offset (s).
%     diag     struct: .Ovec, .scoreCombined, .scoreObject, .offset,
%              .peakcorr (combined corr at offset), .candidates.
%
%   See also LOADFUSISTIMBLOCKS, BUILDFUSISTIMREGRESSORS, PREPROCESSFUSISCAN.

if nargin < 4, opts = struct(); end
opts = setdefault(opts, 'b', []);
opts = setdefault(opts, 'a', []);
opts = setdefault(opts, 'hrfparams', [1.5 10 0.5 1 20 0 16]);
opts = setdefault(opts, 'coarsestep', 1);
opts = setdefault(opts, 'finestep', 0.25);
opts = setdefault(opts, 'finehalfwin', 3);
opts = setdefault(opts, 'searchwin', []);
opts = setdefault(opts, 'prioroffset', NaN);
opts = setdefault(opts, 'priorhalfwin', 90);
opts = setdefault(opts, 'verbose', true);

tframes = double(tframes(:));
nt      = numel(tframes);
RT      = median(diff(tframes));
hrf     = hemodynamicResponse(RT, opts.hrfparams);

on  = blocks.onset(:).';
off = blocks.offset(:).';
selObj = blocks.isObject(:).';

% feasible offsets: keep every block inside the recording window. The
% binding constraints are first-block onset >= scan start and last-block
% offset <= scan end, i.e. O in [t(1)-on(1), t(end)-off(end)] (widened 5 s).
feaslo = tframes(1)   - on(1)    - 5;
feashi = tframes(end) - off(end) + 5;
if ~isempty(opts.searchwin)
    lo = opts.searchwin(1);  hi = opts.searchwin(2);
elseif isfinite(opts.prioroffset)
    % centre on the time0-based prior, clamped to the feasible range
    lo = max(feaslo, opts.prioroffset - opts.priorhalfwin);
    hi = min(feashi, opts.prioroffset + opts.priorhalfwin);
else
    lo = feaslo;  hi = feashi;
end
if lo >= hi, lo = feaslo; hi = feashi; end       % guard a bad prior window

% z-score the data columns ONCE; then corr with a z-scored regressor is a
% single matrix-vector product per candidate offset.
Xz = zscore(X, 0, 1);                 % [nt x nvox]

% ---- coarse scan --------------------------------------------------------
Ovec = lo:opts.coarsestep:hi;
scoreCombined = zeros(size(Ovec));
scoreObject   = zeros(size(Ovec));
for k = 1:numel(Ovec)
    scoreCombined(k) = maxvoxcorr(Xz, buildreg(Ovec(k), true(size(selObj))), nt);
    scoreObject(k)   = maxvoxcorr(Xz, buildreg(Ovec(k), selObj), nt);
end

% ---- disambiguate the 24 s periodicity with the object regressor -------
% candidates = coarse offsets whose combined score is near the global max
peak       = max(scoreCombined);
candidates = Ovec(scoreCombined >= 0.7*peak);
% among them, prefer the one where object blocks also best explain the data
[~, ic]   = max(arrayfun(@(O) maxvoxcorr(Xz, buildreg(O, selObj), nt), candidates));
Ocoarse   = candidates(ic);

% ---- fine refinement around the chosen coarse offset -------------------
Ofine = (Ocoarse-opts.finehalfwin):opts.finestep:(Ocoarse+opts.finehalfwin);
sf    = arrayfun(@(O) maxvoxcorr(Xz, buildreg(O, true(size(selObj))), nt), Ofine);
[peakcorr, jf] = max(sf);
offset = Ofine(jf);

diag = struct('Ovec', Ovec, 'scoreCombined', scoreCombined, ...
    'scoreObject', scoreObject, 'offset', offset, 'peakcorr', peakcorr, ...
    'candidates', candidates, 'prioroffset', opts.prioroffset);

if opts.verbose
    pmsg = '';
    if isfinite(opts.prioroffset)
        pmsg = sprintf(' [prior %+.1f s, delta %.1f s]', ...
            opts.prioroffset, offset - opts.prioroffset);
    end
    fprintf(['alignFusiStimToScan: offset = %+.2f s ', ...
        '(max voxel corr %.3f; searched [%.0f %.0f] s)%s.\n'], ...
        offset, peakcorr, lo, hi, pmsg);
end
if isfinite(opts.prioroffset) && abs(offset - opts.prioroffset) > 60
    warning('alignFusiStimToScan:priorMismatch', ...
        ['Empirical offset %+.1f s differs from the time0 prior %+.1f s ', ...
         'by %.0f s (> 60 s); check clocks for this session.'], ...
        offset, opts.prioroffset, abs(offset - opts.prioroffset));
end

    % --- nested helpers (capture hrf, on, off, tframes, opts) -----------
    function r = buildreg(O, sel)
        box = double(any(tframes >= on(sel)+O & tframes <= off(sel)+O, 2));
        r   = conv(box, hrf); r = r(1:nt);
        if ~isempty(opts.b), r = filtfilt(opts.b, opts.a, r); end
    end
end

% -------------------------------------------------------------------------
function s = maxvoxcorr(Xz, r, nt)
%MAXVOXCORR Largest Pearson correlation between a regressor and any voxel.
%   Xz columns are z-scored (std uses N-1). With r mean-removed and scaled to
%   unit norm, corr(r, col) = (r' * col)/sqrt(N-1), bounded by 1.
r = r - mean(r);
rn = norm(r);
if rn == 0, s = 0; return; end
r = r / rn;
s = max(r.' * Xz) / sqrt(nt-1);
end

% -------------------------------------------------------------------------
function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
