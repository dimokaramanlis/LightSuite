function res = fusiSessionActivationMaps(sesspath, opts)
%FUSISESSIONACTIVATIONMAPS End-to-end activation maps for one fUS session.
%   res = fusiSessionActivationMaps(sesspath) computes Fig.-1C-style visual
%   activation maps (correlation / beta / T-score for combined, object and
%   scrambled blocks, plus an object-vs-scrambled preference map) for a
%   single anesthetized fUS recording, in native session space.
%
%   Pipeline (each step a reusable function):
%     1) load I + time (+ time0) from the *_FUS.mat file (sesspath).
%     2) find and parse the PsychoPy stimulus log   -> loadFusiStimBlocks
%     3) preprocess (dI/I, drift removal, smoothing) -> preprocessFusiScan
%     4) stimulus->scan offset: by default the deterministic key-press trigger
%        (-blocks.triggerPsyTime) is used directly and the hemodynamic lag is
%        left to the HRF; alignFusiStimToScan runs only as a fallback (no key
%        press) or when opts.refinealign is set.
%     5) build HRF-convolved block regressors          -> buildFusiStimRegressors
%     6) fit the voxelwise GLM / correlation            -> computeFusiActivationMaps
%
%   Input:
%     sesspath  full path to a *_FUS.mat file.
%     opts      (optional) struct:
%        .csvpath   stimulus CSV (default: the single *.csv next to sesspath).
%        .finshape  volume shape to reshape I into (default [36 64 54]).
%        .offset    known stimulus->scan offset (s); overrides everything.
%        .refinealign  if true, refine the key-press offset against the data
%                   (fits the ~1 s lag too); default false (key press direct).
%        .hrfparams HRF params (default [1.5 10 0.5 1 20 0 16], the
%                   mapCorrelation values; sampled at the fUS frame rate).
%        .which     which maps to compute (see computeFusiActivationMaps).
%        .timecourses  if true, also return peristimulus block-triggered
%                   responses (fusiPeristimulusTimecourses); default false.
%        .tcwindow  peristimulus window [tpre tpost] s (default [-6 24]).
%        .roithresh combined-corr threshold defining the timecourse ROI
%                   (default 0.2).
%        .preprocess struct forwarded to preprocessFusiScan.
%        .savepath  if set, save res to <savepath>/<session>_actmaps.mat.
%        .verbose   default true.
%
%   Output struct res:
%     .maps       from computeFusiActivationMaps (native-space volumes)
%     .timecourses peristimulus responses (fusiPeristimulusTimecourses) or []
%     .offset     stimulus->scan offset used (s)
%     .blocks     from loadFusiStimBlocks
%     .info       from preprocessFusiScan (fs, filter, volshape, ...)
%     .aligndiag  alignment diagnostics ([] if opts.offset given)
%     .sesspath, .sessionname, .csvpath
%
%   See also LOADFUSISTIMBLOCKS, PREPROCESSFUSISCAN, ALIGNFUSISTIMTOSCAN,
%   BUILDFUSISTIMREGRESSORS, COMPUTEFUSIACTIVATIONMAPS, PLOTFUSIACTIVATIONMONTAGE.

if nargin < 2, opts = struct(); end
opts = setdefault(opts, 'finshape', [36 64 54]);
opts = setdefault(opts, 'hrfparams', [1.5 10 0.5 1 20 0 16]);
opts = setdefault(opts, 'which', {'combined','object','scrambled','preference'});
opts = setdefault(opts, 'preprocess', struct());
opts = setdefault(opts, 'offset', []);
opts = setdefault(opts, 'refinealign', false);
opts = setdefault(opts, 'timecourses', false);   % also return peristimulus TCs
opts = setdefault(opts, 'roithresh', 0.2);        % combined-corr ROI for the TC
opts = setdefault(opts, 'tcwindow', [-6 24]);     % peristimulus window (s)
opts = setdefault(opts, 'savepath', '');
opts = setdefault(opts, 'verbose', true);

[sessdir, sessname] = fileparts(sesspath);
sessname = erase(sessname, '_FUS');

% ---- locate the stimulus CSV ------------------------------------------
if ~isfield(opts, 'csvpath') || isempty(opts.csvpath)
    cand = dir(fullfile(sessdir, '*.csv'));
    if isempty(cand)
        error('fusiSessionActivationMaps:noCSV', ...
            'No stimulus .csv found next to %s', sesspath);
    end
    opts.csvpath = fullfile(cand(1).folder, cand(1).name);
end

if opts.verbose
    fprintf('=== %s ===\n', sessname);
    fprintf('  scan: %s\n  stim: %s\n', sesspath, opts.csvpath);
end

% ---- (1) load the recording -------------------------------------------
S = load(sesspath, 'I', 'time', 'time0');
I = reshape(S.I, [opts.finshape size(S.I, 3)]);
tframes = S.time(:);
time0   = [];
if isfield(S, 'time0'), time0 = S.time0; end   % fUS start wall-clock datenum
clear S;

% ---- (2) stimulus blocks ----------------------------------------------
blocks = loadFusiStimBlocks(opts.csvpath);

% Stimulus->scan offsets from the wall clocks (no fUS signal needed):
%  1) key press: the space press triggers the fUS clock, so the offset is
%     just -triggerPsyTime. The ~1 s hemodynamic lag is left for the HRF to
%     model, so this is used directly (no empirical refinement) by default.
%  2) time0/date: (PsychoPy start - fUS start), good only to ~+/-30 s; used
%     as a coarse prior for the empirical alignment when there is no key press.
keyOffset = NaN; time0Offset = NaN;
if isfinite(blocks.triggerPsyTime), keyOffset = -blocks.triggerPsyTime; end
if ~isempty(time0) && isfinite(blocks.psyStartDatenum)
    time0Offset = (blocks.psyStartDatenum - time0) * 86400;
end
if isfinite(time0Offset)
    prioroffset = time0Offset; priorhalfwin = 90;   % coarse wall-clock prior
else
    prioroffset = NaN;         priorhalfwin = 90;   % full feasible scan
end
if opts.verbose
    fprintf('  offset priors: key-press %+.1f s, time0 %+.1f s\n', ...
        keyOffset, time0Offset);
end

% ---- (3) preprocess ----------------------------------------------------
[X, info] = preprocessFusiScan(I, tframes, opts.preprocess);
clear I;

% ---- (4) determine the stimulus->scan offset --------------------------
% Default: use the deterministic key-press offset directly; the HRF (standard
% mapCorrelation params, sampled at the fUS frame rate) accounts for the
% hemodynamic lag, so no signal-based refinement is done. The empirical
% alignment runs only when there is no key press, or when opts.refinealign is
% set to also fit the ~1 s lag from the data.
aligndiag = [];
if ~isempty(opts.offset)
    offset = opts.offset;
    if opts.verbose, fprintf('  using supplied offset %+.2f s\n', offset); end
elseif isfinite(keyOffset) && ~opts.refinealign
    offset = keyOffset;
    aligndiag = struct('offset', offset, 'keyoffset', keyOffset, ...
        'time0offset', time0Offset, 'method', 'keypress');
    if opts.verbose
        fprintf(['  using key-press offset %+.2f s directly ', ...
            '(HRF models the lag; no refinement).\n'], offset);
    end
    if isfinite(time0Offset) && abs(keyOffset - time0Offset) > 45
        warning('fusiSessionActivationMaps:priorMismatch', ...
            ['key-press offset %+.1f s and time0 estimate %+.1f s differ ', ...
             'by %.0f s; check the clocks for this session.'], ...
            keyOffset, time0Offset, abs(keyOffset - time0Offset));
    end
else
    % no key press (or refinement requested): empirical cross-correlation,
    % centred on the key-press offset when available, else the time0 prior.
    if isfinite(keyOffset), prioroffset = keyOffset; priorhalfwin = 10; end
    aopts = struct('b', info.b, 'a', info.a, 'hrfparams', opts.hrfparams, ...
        'prioroffset', prioroffset, 'priorhalfwin', priorhalfwin, ...
        'verbose', opts.verbose);
    [offset, aligndiag] = alignFusiStimToScan(X, tframes, blocks, aopts);
    if isfinite(keyOffset) && aligndiag.peakcorr < 0.15
        if opts.verbose
            fprintf(['  key-press prior gave weak fit (corr %.2f); ', ...
                'retrying with full scan.\n'], aligndiag.peakcorr);
        end
        aopts.prioroffset = NaN; aopts.priorhalfwin = 90;
        [offset, aligndiag] = alignFusiStimToScan(X, tframes, blocks, aopts);
    end
end

% ---- (5) regressors + (6) maps ----------------------------------------
R    = buildFusiStimRegressors(blocks, tframes, offset, opts.hrfparams);
maps = computeFusiActivationMaps(X, R, info, opts.which);

% ---- (7) peristimulus timecourses (optional) --------------------------
timecourses = [];
if ~isequal(opts.timecourses, false)
    roimask = maps.combined.corr(:) > opts.roithresh;
    timecourses = fusiPeristimulusTimecourses(X, tframes, blocks, offset, ...
        struct('window', opts.tcwindow, 'roimask', roimask));
    if opts.verbose
        fprintf(['  peristimulus timecourses: %d/%d/%d blocks, ', ...
            'ROI %d voxels (combined corr > %.2f)\n'], ...
            timecourses.nblocks, sum(roimask), opts.roithresh);
    end
end

res = struct('maps', maps, 'timecourses', timecourses, 'offset', offset, ...
    'blocks', blocks, 'info', rmfield_safe(info, 'baseline'), ...
    'aligndiag', aligndiag, 'sesspath', string(sesspath), ...
    'sessionname', string(sessname), 'csvpath', string(opts.csvpath));

% ---- optional save -----------------------------------------------------
if ~isempty(opts.savepath)
    makeNewDir(opts.savepath);
    outfile = fullfile(opts.savepath, sprintf('%s_actmaps.mat', sessname));
    save(outfile, '-struct', 'res', '-v7.3');
    if opts.verbose, fprintf('  saved maps to %s\n', outfile); end
end
end

% -------------------------------------------------------------------------
function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function s = rmfield_safe(s, f)
if isfield(s, f), s = rmfield(s, f); end
end
