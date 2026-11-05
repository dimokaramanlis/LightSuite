function out = atlasVoxelGlmMap(mousedata, avsize, sigma, y, mask, opts)
%ATLASVOXELGLMMAP Voxelwise Gaussian GLM of per-animal density maps on a design.
%
%   out = atlasVoxelGlmMap(mousedata, avsize, sigma, y, mask)
%   out = atlasVoxelGlmMap(..., 'Scale', s)
%
%   For every voxel kept by MASK, fits across animals the Gaussian
%   (identity-link) GLM
%
%       map_a(v) = b0(v) + b_y(v)' * y_a + b_s(v) * s_a + e_a(v),
%       e_a(v) ~ N(0, sigma2(v))
%
%   where map_a is the smoothed map of animal a. Smoothed counts are sums of many
%   points, so a Gaussian error model is reasonable. The per-animal scale (e.g.
%   total cell yield) enters as its own regressor rather than rescaling the
%   data, so it gets its own coefficient and t-map.
%
%   With 'Transform','log' the response is log(map + c) and the scale enters
%   as log(Scale), i.e. a multiplicative model of the density:
%
%       log(map_a(v) + c) = b0(v) + b_y(v)' * y_a + b_s(v) * log(s_a) + e_a(v)
%
%   This matches TESTBRAINAREACOVARIATELM at the voxel level. b_y is then a
%   log fold change per unit of y. With a single regressor of interest, its
%   t-test is the exact (F) version of the likelihood-ratio test of adding y.
%
%   The design matrix is the same at every voxel, so the fit is closed form.
%   With a thin QR, X = Q*R, the only per-voxel quantities needed are Q'*map
%   (p numbers) and sum(map.^2). They are accumulated animal by animal while
%   streaming z-slabs, exactly as in ATLASVOXELCORRMAP, so speed and memory
%   are about the same as the correlation map.
%
%   INPUTS (same conventions as ATLASVOXELCORRMAP)
%     mousedata : cell array, one entry per animal: [Npoints x 3] atlas-space
%                 points ([x y z], used as (:,[2 1 3])) or a function handle
%                 @(z0,z1) returning [avsize(1) avsize(2) z1-z0+1] on the
%                 OUTPUT z-grid.
%     avsize    : [Ny Nx Nz] atlas volume size.
%     sigma     : Gaussian sigma in voxels, scalar or [sy sx sz].
%     y         : [Nmice x K] regressor(s) of interest. All K columns enter the
%                 same model jointly (multiple regression).
%     mask      : logical array of the OUTPUT volume size, true = keep.
%                 [] keeps every voxel.
%
%   NAME-VALUE OPTIONS
%     'Scale'            []   - [Nmice x 1] per-animal scale (e.g.
%                               cellfun(@(c) size(c,1), mousedata)), added as a
%                               regressor as is. [] = no scale term.
%     'Transform'        'none' - 'none' fits the smoothed map as is; 'log'
%                               fits log(map + c) and uses log(Scale).
%     'Pseudocount'      1    - c for 'log', in effective cells under the
%                               kernel. A smoothed map is a density per voxel,
%                               so c is converted to map units as
%                               Pseudocount / ((4*pi)^1.5 * prod(sigma)), the
%                               Poisson-equivalent cell count of the kernel.
%                               Assumes maps in cells per voxel (point input,
%                               or loaders on the same scale).
%     'CenterRegressors' true - centre y and Scale so the intercept is the
%                               map of the average animal. Slopes and their t
%                               are unaffected.
%     'Symmetrize'       true - mirror the hemispheres; output is [Ny Nx Nz/2].
%     'SlabThickness'    0    - slab depth in voxels; 0 = pick from memory budget.
%     'MemoryBudgetGB'   2
%     'UseGPU'           auto
%     'Tail'             'both' - 'both' | 'right' | 'left' (all coefficients).
%     'OutputFormat'     'masked' - 'full' also returns NaN-filled volumes
%                               betaVol / tVol / pVol for the first y column.
%     'Verbose'          true
%
%   OUTPUT (struct)
%     out.idx        [Nvox x 1] uint32 linear indices into the output volume
%     out.coefNames  {1 x p} 'Intercept', 'y1'..'yK', 'Scale'
%     out.beta       [Nvox x p] single, coefficients (columns as coefNames)
%     out.t          [Nvox x p] single, t-statistics
%     out.p          [Nvox x p] single, p-values, df = Nmice - p
%                    (NaN for the intercept)
%     out.sigma2     [Nvox x 1] single, residual variance (RSS/df)
%     out.X          [Nmice x p] design matrix that was fitted
%     out.volumeSize, out.n, out.df, out.tail, out.sigma
%     out.transform, out.pseudocount (in map units, 0 for 'none')
%     out.betaVol / out.tVol / out.pVol when 'OutputFormat' is 'full'
%
%   EXAMPLE
%     ncells = cellfun(@(c) size(c,1), mousecells);
%     out    = atlasVoxelGlmMap(mousecells, avsize, 4, behaviorScore, mask, ...
%                  'Scale', ncells, 'UseGPU', true);
%     iy     = strcmp(out.coefNames, 'y1');
%     tvol   = nan(out.volumeSize, 'single');
%     tvol(out.idx) = out.t(:, iy);
%
%   See also ATLASVOXELCORRMAP, ATLASCELLMAP, FITLM.

arguments
    mousedata cell
    avsize (1,3) double {mustBePositive}
    sigma (1,:) double {mustBeNonnegative}
    y double
    mask = []
    opts.Scale double = []
    opts.Transform (1,:) char {mustBeMember(opts.Transform, {'none','log'})} = 'none'
    opts.Pseudocount (1,1) double {mustBePositive} = 1
    opts.CenterRegressors (1,1) logical = true
    opts.Symmetrize (1,1) logical = true
    opts.SlabThickness (1,1) double = 0
    opts.MemoryBudgetGB (1,1) double = 2
    opts.UseGPU (1,1) logical = hasGPU()
    opts.Tail (1,:) char {mustBeMember(opts.Tail, {'both','right','left'})} = 'both'
    opts.OutputFormat (1,:) char {mustBeMember(opts.OutputFormat, {'masked','full'})} = 'masked'
    opts.Verbose (1,1) logical = true
end

%==========================================================================
% 1. Geometry and design matrix
%==========================================================================
nA = numel(mousedata);
if size(y, 1) ~= nA
    error('atlasVoxelGlmMap:ySize', 'y must have %d rows (one per animal).', nA);
end
K = size(y, 2);

outsize = avsize;
if opts.Symmetrize
    if mod(avsize(3), 2) ~= 0
        error('atlasVoxelGlmMap:oddZ', 'Symmetrize requires an even avsize(3).');
    end
    outsize(3) = avsize(3)/2;
end
Ny = outsize(1); Nx = outsize(2); Nz = outsize(3);

Xc        = double(y);
coefNames = [{'Intercept'}, arrayfun(@(k) sprintf('y%d', k), 1:K, 'UniformOutput', false)];
if ~isempty(opts.Scale)
    sc = double(opts.Scale(:));
    assert(numel(sc) == nA, 'Scale must have one entry per animal.');
    if strcmp(opts.Transform, 'log')
        assert(all(sc > 0), 'Scale must be positive for the log transform.');
        sc = log(sc);
    end
    Xc        = [Xc, sc];
    coefNames = [coefNames, {'Scale'}];
end
if opts.CenterRegressors
    Xc = Xc - mean(Xc, 1);
end
X  = [ones(nA, 1), Xc];
p  = size(X, 2);
df = nA - p;
if df < 1
    error('atlasVoxelGlmMap:tooFewAnimals', ...
        '%d animals for %d coefficients: need at least %d.', nA, p, p+1);
end

% Thin QR: beta = R \ (Q'v), RSS = v'v - ||Q'v||^2, cov(beta) = sigma2*inv(X'X)
[Q, R] = qr(X, 0);
if rank(R) < p
    error('atlasVoxelGlmMap:rankDeficient', 'Design matrix is rank deficient.');
end
Rinv  = R \ eye(p);
cdiag = sum(Rinv.^2, 2)';                   % diag(inv(X'X)), 1 x p
Qs    = single(Q);
RinvT = single(Rinv');

if isscalar(sigma), sigma = repmat(sigma, 1, 3); end
zrad = ceil(2*sigma(3));                    % imgaussfilt3 default FilterSize radius

% log transform: pseudocount in map units. A Gaussian-smoothed Poisson
% density has the variance of a count of (4*pi)^1.5*prod(sigma)*density
% cells, so this is "Pseudocount cells" on that effective-count scale.
dolog = strcmp(opts.Transform, 'log');
if dolog
    keff = (4*pi)^1.5 * prod(max(sigma, 0.5));
    pc   = single(opts.Pseudocount / keff);
else
    pc   = single(0);
end

%==========================================================================
% 2. Mask and per-animal point prep
%==========================================================================
haveMask = ~isempty(mask);
if haveMask
    assert(isequal(size(mask), outsize), ...
        'mask must be %dx%dx%d (the output volume size).', Ny, Nx, Nz);
    mask = logical(mask);
end

pts = cell(nA, 1);
for a = 1:nA
    if isa(mousedata{a}, 'function_handle')
        pts{a} = mousedata{a};
    else
        pts{a} = prepPoints(mousedata{a}, outsize, opts.Symmetrize);
    end
end

%==========================================================================
% 3. Slab thickness and GPU
%==========================================================================
if opts.SlabThickness > 0
    nzSlab = round(opts.SlabThickness);
else
    % per slab: working volume (+ smoothing copy, buffer), plus p+2 masked accumulators
    bytesPerSlice = Ny*Nx*4;
    nArrays       = p + 4;
    nzSlab        = floor(opts.MemoryBudgetGB*2^30/bytesPerSlice/nArrays) - 2*zrad;
end
nzSlab = max(1, min(Nz, nzSlab));

useGPU = opts.UseGPU;
if useGPU && ~hasGPU()
    warning('atlasVoxelGlmMap:noGPU', 'No GPU available, running on the CPU.');
    useGPU = false;
end
if useGPU
    Qs = gpuArray(Qs); RinvT = gpuArray(RinvT);
end

%==========================================================================
% 4. Slab loop
%==========================================================================
z0list = 1:nzSlab:Nz;
nS     = numel(z0list);
idxC   = cell(nS, 1);
bC     = cell(nS, 1);
rssC   = cell(nS, 1);
tstart = tic;

for s = 1:nS
    z0    = z0list(s);
    z1    = min(Nz, z0 + nzSlab - 1);
    zl    = max(1,  z0 - zrad);
    zh    = min(Nz, z1 + zrad);
    nzb   = zh - zl + 1;
    nzOut = z1 - z0 + 1;
    keepZ = (z0 - zl + 1) + (0:nzOut-1);

    if haveMask
        sel = mask(:,:,z0:z1);
        if ~any(sel(:)), continue; end
    else
        sel = true(Ny, Nx, nzOut);
    end
    if useGPU, selg = gpuArray(sel); else, selg = sel; end
    nvS   = nnz(sel);
    allIn = nvS == numel(sel);

    % Sums are taken on d_a = map_a - map_1. With an intercept in the model
    % this leaves residuals and slopes unchanged and only shifts the
    % intercept by map_1, but it avoids single-precision cancellation in
    % sum(v.^2) - ||Q'v||^2 when maps sit on a large baseline.
    QtD = zeros(nvS, p, 'single');
    SDD = zeros(nvS, 1, 'single');
    if useGPU, QtD = gpuArray(QtD); SDD = gpuArray(SDD); end

    for a = 1:nA
        if isa(pts{a}, 'function_handle')
            V = single(pts{a}(zl, zh));
            assert(isequal(size(V), [Ny Nx nzb]), ...
                'Loader for animal %d returned %s, expected [%d %d %d].', ...
                a, mat2str(size(V)), Ny, Nx, nzb);
        else
            V = buildSlab(pts{a}, Ny, Nx, zl, zh);
        end
        if useGPU, V = gpuArray(V); end

        V = imgaussfilt3(V, sigma, 'Padding', 'replicate');
        V = V(:,:,keepZ);
        if allIn
            v = V(:);                       % no copy when every voxel is kept
        else
            v = V(selg);                    % only in-mask voxels from here on
        end
        if dolog
            v = log(max(v, 0) + pc);        % max: smoothing can leave tiny negatives
        end

        if a == 1
            ref = v;
            continue                        % d_1 = 0 contributes nothing
        end
        d   = v - ref;
        for k = 1:p                         % column-wise, updated in place
            QtD(:,k) = QtD(:,k) + d * Qs(a,k);
        end
        SDD = SDD + d.*d;
    end

    bS      = QtD * RinvT;                  % beta of d, nvS x p
    bS(:,1) = bS(:,1) + ref;                % undo the shift on the intercept
    rssS    = max(SDD - sum(QtD.^2, 2), 0);

    idxC{s} = uint32(find(sel)) + uint32(z0-1)*uint32(Ny*Nx);
    bC{s}   = gather(bS);
    rssC{s} = gather(rssS);

    if opts.Verbose
        fprintf('atlasVoxelGlmMap: slab %d/%d (z %d-%d), %d voxels, %.1f s\n', ...
            s, nS, z0, z1, nvS, toc(tstart));
    end
    clear QtD SDD ref d v V bS rssS selg
end

%==========================================================================
% 5. Assemble; t and p only over the kept voxels
%==========================================================================
out.idx       = cat(1, idxC{:});
out.coefNames = coefNames;
beta          = double(cat(1, bC{:}));
sigma2        = double(cat(1, rssC{:})) / df;

se = sqrt(sigma2 .* cdiag);                 % Nvox x p
t  = beta ./ se;
t(sigma2 <= 0, :) = NaN;                    % no residual variance: undefined

% p in single (and on the GPU if used): this is the slowest step otherwise.
% The intercept is skipped, a p-value for "map == 0" is not meaningful.
pv = nan(size(t), 'single');
tt = single(t(:, 2:end));
if useGPU, tt = gpuArray(tt); end
% NaN t (no residual variance, e.g. empty voxels after the log) makes the
% GPU betainc error out, so feed it a valid x there and NaN the p afterwards
xb    = df ./ (df + tt.^2);
xb(~isfinite(xb)) = 1;
xb    = min(max(xb, 0), 1);
pboth = betainc(xb, single(df/2), single(0.5));   % == 2*tcdf(-|t|, df)
clear xb
switch opts.Tail
    case 'both'
        pp = pboth;
    case 'right'
        pp = pboth/2;  pp(tt <= 0) = 1 - pboth(tt <= 0)/2;
    case 'left'
        pp = pboth/2;  pp(tt >= 0) = 1 - pboth(tt >= 0)/2;
end
pp(~isfinite(tt)) = NaN;
pv(:, 2:end) = gather(pp);
clear tt pboth pp

out.beta       = single(beta);
out.t          = single(t);
out.p          = pv;
out.sigma2     = single(sigma2);
out.X          = X;
out.n          = nA;
out.df         = df;
out.tail       = opts.Tail;
out.sigma      = sigma;
out.volumeSize  = outsize;
out.transform   = opts.Transform;
out.pseudocount = double(pc);

if strcmp(opts.OutputFormat, 'full')
    out.betaVol = nan(outsize, 'single');
    out.tVol    = nan(outsize, 'single');
    out.pVol    = nan(outsize, 'single');
    out.betaVol(out.idx) = out.beta(:,2);
    out.tVol(out.idx)    = out.t(:,2);
    out.pVol(out.idx)    = out.p(:,2);
    if K > 1
        warning('atlasVoxelGlmMap:fullFirstColumn', ...
            'Full volumes were filled with the first column of y only.');
    end
end

if opts.Verbose
    fprintf('atlasVoxelGlmMap: %d voxels, %d animals, %d coefs, done in %.1f s\n', ...
        numel(out.idx), nA, p, toc(tstart));
end
end

%==========================================================================
function P = prepPoints(cellcoords, outsize, symflag)
% Convert a point list into slice-sorted subscripts + weights, once per animal.

s = round(double(cellcoords(:, [2 1 3])));     % [row col page], as atlasCellMap
w = ones(size(s,1), 1, 'single');

if symflag
    Nmid    = outsize(3);
    hi      = s(:,3) > Nmid;
    s(hi,3) = 2*Nmid - s(hi,3) + 1;            % mirror onto the kept half
    w       = w/2;                             % atlasCellMap averages the halves
end

ok = all(s >= 1, 2) & s(:,1) <= outsize(1) & s(:,2) <= outsize(2) & s(:,3) <= outsize(3);
s  = s(ok,:);
w  = w(ok);

[z, ord] = sort(s(:,3));
P.rc   = uint32(s(ord,1) + (s(ord,2)-1)*outsize(1));
P.z    = uint32(z);
P.w    = w(ord);
P.edge = [0; cumsum(accumarray(z, 1, [outsize(3) 1]))];
end

%==========================================================================
function V = buildSlab(P, Ny, Nx, zl, zh)
% Rasterize one animal's points into a [Ny Nx zh-zl+1] single volume.

nzb = zh - zl + 1;
i0  = P.edge(zl) + 1;
i1  = P.edge(zh+1);
if i1 < i0
    V = zeros(Ny, Nx, nzb, 'single');
    return
end
sel = i0:i1;
lin = double(P.rc(sel)) + (double(P.z(sel)) - zl)*(Ny*Nx);
V   = accumarray(lin, double(P.w(sel)), [Ny*Nx*nzb, 1]);
V   = reshape(single(V), Ny, Nx, nzb);
end

%==========================================================================
function tf = hasGPU()
try
    tf = gpuDeviceCount > 0;
catch
    tf = false;
end
end
