function blocks = loadFusiStimBlocks(csvpath)
%LOADFUSISTIMBLOCKS Parse a PsychoPy blockedTrials CSV into stimulus blocks.
%   blocks = loadFusiStimBlocks(csvpath) reads one *_blockedTrials_*.csv log
%   from the anesthetized fUS visual-object experiment and returns the block
%   structure of the run. The run is a block design: 12 s gray followed by
%   12 s of images (12 images, ~1 s each), repeated for 64 blocks (32 object,
%   32 texture-scrambled) in pseudorandom order (see Fig. 1B of the paper).
%
%   Each image row logs its onset in the PsychoPy clock in column
%   "image.started"; the block a row belongs to is identified by the block
%   file name in "condsFile" (names containing "texture" are the scrambled
%   controls). Consecutive image onsets within a block are ~1 s apart, so
%   blocks are split wherever the gap between successive onsets exceeds 3 s.
%
%   IMPORTANT: the returned times are in the PsychoPy clock, which does NOT
%   share an origin with the fUS acquisition clock (FUS.mat "time"). The
%   per-session offset between the two must be estimated separately with
%   alignFusiStimToScan before the blocks can be placed on the scan grid.
%
%   Output struct blocks:
%     .onset      [nBlk x 1] first image onset of each block (PsychoPy s)
%     .offset     [nBlk x 1] last image onset + 0.5 s (block stim end, s)
%     .isObject   [nBlk x 1] logical, true for object blocks
%     .isScrambled[nBlk x 1] logical, true for texture-scrambled blocks
%     .category   [nBlk x 1] string, block family (e.g. "Birds","Mice",...)
%     .condsFile  [nBlk x 1] string, raw condsFile of each block
%     .nImages    [nBlk x 1] number of image onsets in the block
%     .psyStartDatenum  PsychoPy session start as a MATLAB datenum, parsed
%                    from the CSV "date" field (minute precision), or NaN.
%                    With FUS.mat "time0" (fUS start datenum) this gives a
%                    coarse (~+/-30 s) prior on the stimulus->scan offset.
%     .triggerPsyTime  absolute time (PsychoPy clock, s) of the run-start key
%                    press = key_resp.started + key_resp.rt, or NaN. This is
%                    the hardware trigger that zeroes the fUS clock, so the
%                    stimulus->scan offset is -triggerPsyTime to ~1 s (the
%                    residual being hemodynamic lag). Much tighter than the
%                    time0/date prior; used by alignFusiStimToScan.
%     .csvpath    source file
%
%   See also ALIGNFUSISTIMTOSCAN, BUILDFUSISTIMREGRESSORS.

T = readtable(csvpath, 'VariableNamingRule', 'preserve');

% image onset column and the block-identity column
ist  = T.('image.started');
cond = string(T.('condsFile'));

% keep only real image-presentation rows (numeric onset, non-empty block)
if ~isnumeric(ist), ist = str2double(string(ist)); end
valid = ~isnan(ist) & ~ismissing(cond) & cond ~= "";
ist   = ist(valid);
cond  = cond(valid);

% sort by time so the gap-based block segmentation is well defined
[ist, ord] = sort(ist);
cond       = cond(ord);

if isempty(ist)
    error('loadFusiStimBlocks:noImages', ...
        'No image.started rows found in %s', csvpath);
end

% segment into blocks: a new block starts wherever the onset gap jumps (>3 s)
newblk = [true; diff(ist) > 3];
bstart = find(newblk);
bend   = [bstart(2:end) - 1; numel(ist)];

condFirst = cond(bstart);
isScr     = contains(condFirst, 'texture', 'IgnoreCase', true);

blocks = struct();
blocks.onset       = ist(bstart);
blocks.offset      = ist(bend) + 0.5;                 % + last 0.5 s gray tail
blocks.isScrambled = isScr;
blocks.isObject    = ~isScr;
blocks.category    = extractBefore(erase(condFirst, '.csv'), 'Block');
blocks.condsFile   = condFirst;
blocks.nImages     = bend - bstart + 1;
blocks.psyStartDatenum = parsePsyStart(T);
blocks.triggerPsyTime  = parseTrigger(T);
blocks.csvpath     = string(csvpath);

fprintf(['loadFusiStimBlocks: %d blocks (%d object, %d scrambled) ', ...
    'spanning %.1f-%.1f s (PsychoPy clock).\n'], numel(blocks.onset), ...
    sum(blocks.isObject), sum(blocks.isScrambled), ...
    blocks.onset(1), blocks.offset(end));
end

% -------------------------------------------------------------------------
function dn = parsePsyStart(T)
%PARSEPSYSTART PsychoPy session start (datenum) from the CSV "date" field.
%   The field looks like "2022_Feb_21_0932" (minute precision). Returns NaN
%   if the column is absent or unparseable.
dn = NaN;
if ~any(strcmp('date', T.Properties.VariableNames)), return; end
ds = string(T.('date'));
ds = ds(find(~ismissing(ds) & ds ~= "", 1));
if isempty(ds), return; end
try
    dn = datenum(datetime(ds, 'InputFormat', 'yyyy_MMM_dd_HHmm'));
catch
    dn = NaN;
end
end

% -------------------------------------------------------------------------
function tp = parseTrigger(T)
%PARSETRIGGER Absolute PsychoPy-clock time (s) of the run-start key press.
%   Returns key_resp.started + key_resp.rt of the first logged press, or NaN.
tp = NaN;
vn = T.Properties.VariableNames;
if ~all(ismember({'key_resp.started', 'key_resp.rt'}, vn)), return; end
st = tonum(T.('key_resp.started'));
rt = tonum(T.('key_resp.rt'));
k  = find(isfinite(st) & isfinite(rt), 1);
if ~isempty(k), tp = st(k) + rt(k); end
end

% -------------------------------------------------------------------------
function x = tonum(x)
if ~isnumeric(x), x = str2double(string(x)); end
x = x(:);
end
