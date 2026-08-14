function R = buildFusiStimRegressors(blocks, tframes, offset, hrfparams)
%BUILDFUSISTIMREGRESSORS Block-design fUS regressors on the scan time grid.
%   R = buildFusiStimRegressors(blocks, tframes, offset) turns the stimulus
%   blocks from loadFusiStimBlocks into predicted fUS response regressors,
%   sampled at the frame times tframes (FUS.mat "time", seconds), for the
%   three contrasts used in Fig. 1:
%       combined   - every image block (object + scrambled), Fig. 1C
%       object     - object blocks only
%       scrambled  - texture-scrambled blocks only
%
%   Each regressor is a 0/1 block boxcar (1 while a block's images are on the
%   screen) convolved with a single hemodynamic response function, so it is
%   the predicted power-Doppler time course of a voxel that follows that
%   stimulus. The boxcars are built on the (possibly non-uniform) tframes
%   grid directly, so no resampling of the data is required.
%
%   offset (seconds, default 0) is the stimulus->scan clock offset from
%   alignFusiStimToScan: a block written in the PsychoPy clock at time p is
%   placed on the scan grid at t = p + offset. It MUST be supplied (or the
%   maps will be aligned to noise) because the PsychoPy and fUS clocks do not
%   share an origin.
%
%   hrfparams (default [1.5 10 0.5 1 20 0 16], the lab's fUS HRF used in
%   mapCorrelation) are passed straight to hemodynamicResponse, sampled at
%   RT = median(diff(tframes)).
%
%   Output struct R:
%     .combined/.object/.scrambled  [nt x 1] HRF-convolved regressors (raw,
%                                    not z-scored; zero-mean removed)
%     .boxCombined/.boxObject/.boxScrambled  [nt x 1] the raw 0/1 boxcars
%     .hrf        the HRF kernel used
%     .offset     the offset applied
%     .tframes    the frame times used
%
%   See also LOADFUSISTIMBLOCKS, ALIGNFUSISTIMTOSCAN, HEMODYNAMICRESPONSE,
%   COMPUTEFUSIACTIVATIONMAPS.

if nargin < 3 || isempty(offset),    offset = 0; end
if nargin < 4 || isempty(hrfparams), hrfparams = [1.5 10 0.5 1 20 0 16]; end

tframes = double(tframes(:));
RT      = median(diff(tframes));
hrf     = hemodynamicResponse(RT, hrfparams);

on  = blocks.onset(:).'  + offset;   % shift blocks into the scan clock
off = blocks.offset(:).' + offset;

% boxcar: frame is "on" if its time falls inside any selected block window
boxcar = @(sel) double(any(tframes >= on(sel) & tframes <= off(sel), 2));
% convolve with the HRF and crop back to nt (causal)
convhrf = @(b) local_conv(b, hrf);

selObj = blocks.isObject(:).';
selScr = blocks.isScrambled(:).';
selAll = true(size(selObj));

R = struct();
R.boxCombined  = boxcar(selAll);
R.boxObject    = boxcar(selObj);
R.boxScrambled = boxcar(selScr);
R.combined     = convhrf(R.boxCombined);
R.object       = convhrf(R.boxObject);
R.scrambled    = convhrf(R.boxScrambled);
R.hrf          = hrf;
R.offset       = offset;
R.tframes      = tframes;
end

% -------------------------------------------------------------------------
function y = local_conv(x, h)
%LOCAL_CONV Causal convolution cropped to the length of x, mean removed.
y = conv(x, h);
y = y(1:numel(x));
y = y - mean(y);
end
