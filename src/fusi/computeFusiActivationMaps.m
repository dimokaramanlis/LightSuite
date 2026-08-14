function maps = computeFusiActivationMaps(X, R, info, which)
%COMPUTEFUSIACTIVATIONMAPS Voxelwise fUS activation maps (corr / beta / T).
%   maps = computeFusiActivationMaps(X, R, info) fits the block-design
%   regressors R (from buildFusiStimRegressors) to the preprocessed data X
%   (from preprocessFusiScan) and returns, for each voxel, the maps used in
%   Fig. 1:
%     - Pearson correlation with the stimulus regressor,
%     - GLM slope (beta), and
%     - GLM T-score (beta / standard error),
%   for three single-regressor models:
%     .combined   all image blocks (object + scrambled) -> Fig. 1C T-scores
%     .object     object blocks only
%     .scrambled  texture-scrambled blocks only
%   and, from one JOINT model with object and scrambled as separate columns,
%   an object-vs-scrambled preference map (voxelwise analogue of VOSIT):
%     .preference.objBeta / .scrBeta / .objMinusScr (+ _t / _p)
%
%   The GLM is solved in a single vectorized ordinary-least-squares pass
%   (one design shared by all voxels), which is identical to calling glmfit
%   per voxel with a normal distribution but ~1000x faster.
%
%   Input:
%     X     [nt x nvox] preprocessed data (time down columns).
%     R     regressor struct from buildFusiStimRegressors (fields
%           .combined/.object/.scrambled, raw HRF-convolved).
%     info  struct from preprocessFusiScan, used for .b/.a (to filter the
%           regressors exactly like the data) and .volshape (to reshape maps).
%     which (optional) cellstr subset of {'combined','object','scrambled',
%           'preference'} to compute (default all).
%
%   Output struct maps:
%     .combined/.object/.scrambled : each with fields corr, beta, tscore,
%                                    pval (volumes if info.volshape set, else
%                                    [nvox x 1]).
%     .preference : objBeta, scrBeta, objMinusScr, objMinusScr_t,
%                   objMinusScr_p.
%     .df, .volshape, .offset (from R).
%
%   NOTE on statistics: t/p are single-session, uncorrected, and do not model
%   temporal autocorrelation (smoothing inflates effective df), so p-values
%   are anti-conservative. Use them to threshold/visualize a session; formal
%   significance (Fig. 1) needs the cross-session mixed-effects + FDR model.
%
%   See also PREPROCESSFUSISCAN, BUILDFUSISTIMREGRESSORS, MAPCORRELATION.

if nargin < 4 || isempty(which)
    which = {'combined','object','scrambled','preference'};
end
nt       = size(X, 1);
volshape = [];
if isfield(info, 'volshape'), volshape = info.volshape; end
haveHP   = isfield(info,'b') && ~isempty(info.b);

% regressors filtered exactly like the data, then z-scored
prep = @(r) zscore(filtreg(r, info, haveHP));
rComb = prep(R.combined);
rObj  = prep(R.object);
rScr  = prep(R.scrambled);

maps = struct();
maps.volshape = volshape;
maps.offset   = R.offset;

% ---- single-regressor models (corr + beta + T) ------------------------
if any(strcmp(which, 'combined'))
    maps.combined  = singleMap(X, rComb, nt, volshape);
end
if any(strcmp(which, 'object'))
    maps.object    = singleMap(X, rObj, nt, volshape);
end
if any(strcmp(which, 'scrambled'))
    maps.scrambled = singleMap(X, rScr, nt, volshape);
end

% ---- joint object+scrambled model -> preference contrast --------------
if any(strcmp(which, 'preference'))
    D = [ones(nt,1) rObj rScr];
    C = [0 0 0; 1 0 1; 0 1 -1];           % cols: objBeta, scrBeta, obj-scr
    [est, t, p, df] = glmContrasts(X, D, C);
    pref = struct();
    pref.objBeta       = rs(est(1,:), volshape);
    pref.scrBeta       = rs(est(2,:), volshape);
    pref.objMinusScr   = rs(est(3,:), volshape);
    pref.objMinusScr_t = rs(t(3,:),   volshape);
    pref.objMinusScr_p = rs(p(3,:),   volshape);
    maps.preference    = pref;
    maps.df            = df;
else
    maps.df = nt - 2;
end
end

% -------------------------------------------------------------------------
function m = singleMap(X, r, nt, volshape)
%SINGLEMAP corr, beta, T and p for a one-regressor (+intercept) model.
c = corr(X, r);                              % Pearson corr per voxel [nvox x 1]
D = [ones(nt,1) r];
[est, t, p, df] = glmContrasts(X, D, [0;1]); %#ok<ASGLU>
m = struct('corr', rs(c, volshape), 'beta', rs(est(1,:), volshape), ...
    'tscore', rs(t(1,:), volshape), 'pval', rs(p(1,:), volshape));
end

% -------------------------------------------------------------------------
function [est, t, p, df] = glmContrasts(X, D, C)
%GLMCONTRASTS Vectorized OLS: contrasts C'*beta with T-scores over all voxels.
%   X [nt x nvox], D [nt x k], C [k x nc]. Returns est/t/p [nc x nvox].
[nt, k] = size(D);
df      = nt - k;
iDtD    = inv(D.' * D);                 %#ok<MINV> small k, reused below
B       = iDtD * (D.' * X);             % [k x nvox]
resid   = X - D * B;
sigma2  = sum(resid.^2, 1) / df;        % [1 x nvox]
est     = C.' * B;                      % [nc x nvox]
cvar    = diag(C.' * iDtD * C);         % [nc x 1] design part of the variance
se      = sqrt(cvar * sigma2);          % [nc x nvox]
t       = est ./ se;
p       = 2 * tcdf(-abs(t), df);
end

% -------------------------------------------------------------------------
function r = filtreg(r, info, haveHP)
%FILTREG Apply the data's high-pass to a regressor so design matches data.
r = r(:);
if haveHP, r = filtfilt(info.b, info.a, r); end
end

% -------------------------------------------------------------------------
function v = rs(x, volshape)
%RS Reshape a [1 x nvox] or [nvox x 1] map to a volume if volshape is known.
x = x(:);
if isempty(volshape), v = x; else, v = reshape(x, volshape); end
end
