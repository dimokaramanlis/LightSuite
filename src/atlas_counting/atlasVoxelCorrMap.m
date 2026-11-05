function out = atlasVoxelCorrMap(mousedata, avsize, sigma, y, mask, opts)
%ATLASVOXELCORRMAP Voxelwise correlation of per-animal density maps with a vector.
%
%   out = atlasVoxelCorrMap(mousedata, avsize, sigma, y, mask)
%
%   For every voxel kept by MASK, computes the Pearson correlation across animals
%   between the smoothed per-animal map and the Nmice-element vector Y, plus a
%   p-value. The volume is never held in memory for all animals at once: it is
%   processed in z-slabs with a buffer of ceil(2*sigma) slices on each side, so
%   the smoothing is identical to running IMGAUSSFILT3 on the full volume.
%
%   INPUTS (same conventions as ATLASCELLMAP)
%     mousedata : cell array, one entry per animal. Each entry is either
%                 - an [Npoints x 3] array of cell coordinates in atlas space,
%                   columns ordered [x y z]; they are used as (:,[2 1 3]),
%                   exactly like atlasCellMap, or
%                 - a function handle @(z0,z1) returning an array of size
%                   [avsize(1) avsize(2) z1-z0+1] (use this for intensity volumes
%                   or anything already on disk). Slabs are requested on the
%                   OUTPUT z-grid, i.e. after symmetrization.
%     avsize    : [Ny Nx Nz] atlas volume size, as passed to atlasCellMap.
%     sigma     : Gaussian sigma in voxels, scalar or [sy sx sz].
%     y         : [Nmice x 1] vector to correlate against. May be [Nmice x K] to
%                 test K vectors at once (r and p come back with K columns).
%     mask      : logical array of the OUTPUT volume size (see out.volumeSize),
%                 true = keep. [] keeps every voxel.
%
%   NAME-VALUE OPTIONS
%     'Symmetrize'     true  - mirror the two hemispheres like atlasCellMap; the
%                              output volume is then [Ny Nx Nz/2].
%     'Scale'          []    - [Nmice x 1] per-animal multiplicative factor
%                              applied before smoothing, e.g.
%                              1./cellfun(@(c) size(c,1), mousedata) to remove
%                              differences in overall cell yield. This is NOT the
%                              averaging weight of atlasCellMap: a per-animal
%                              scale genuinely changes the correlations.
%     'SlabThickness'  0     - slab depth in voxels; 0 = pick from memory budget.
%     'MemoryBudgetGB' 2     - budget used to pick the slab thickness.
%     'UseGPU'         auto  - run smoothing and accumulation on the GPU.
%     'Tail'           'both'- 'both' | 'right' | 'left'.
%     'OutputFormat'   'masked' - 'masked' returns vectors over in-mask voxels
%                              plus their linear indices; 'full' additionally
%                              expands r and p into NaN-filled volumes (4 bytes
%                              per voxel per map, i.e. 2 GB each at 1000x1000x500).
%     'Verbose'        true
%
%   OUTPUT (struct)
%     out.idx        [Nvox x 1] uint32 linear indices into the output volume
%     out.r          [Nvox x K] single, correlation
%     out.p          [Nvox x K] single, p-value, df = Nmice-2
%     out.meanMap    [Nvox x 1] single, across-animal mean of the smoothed maps
%     out.stdMap     [Nvox x 1] single, across-animal std (normalized by N-1)
%     out.volumeSize size of the output volume
%     out.n, out.df, out.tail, out.sigma
%     out.rVol / out.pVol when 'OutputFormat' is 'full'
%
%   EXAMPLE
%     mask = annotationVolume > 0;                      % brain voxels only
%     mask = mask(:,:,1:size(mask,3)/2);                % symmetrized half
%     nrm  = 1./cellfun(@(c) size(c,1), mousecells);
%     out  = atlasVoxelCorrMap(mousecells, avsize, 4, behaviorScore, mask, ...
%                'Scale', nrm, 'UseGPU', true);
%     sig  = out.p < 0.05;
%     rvol = nan(out.volumeSize, 'single');
%     rvol(out.idx(sig)) = out.r(sig);
%
%   See also ATLASCELLMAP, ATLASCELLMAPVARIANCE, IMGAUSSFILT3.

arguments
    mousedata cell
    avsize (1,3) double {mustBePositive}
    sigma (1,:) double {mustBeNonnegative}
    y double
    mask = []
    opts.Symmetrize (1,1) logical = true
    opts.Scale double = []
    opts.SlabThickness (1,1) double = 0
    opts.MemoryBudgetGB (1,1) double = 2
    opts.UseGPU (1,1) logical = hasGPU()
    opts.Tail (1,:) char {mustBeMember(opts.Tail, {'both','right','left'})} = 'both'
    opts.OutputFormat (1,:) char {mustBeMember(opts.OutputFormat, {'masked','full'})} = 'masked'
    opts.Verbose (1,1) logical = true
end

%==========================================================================
% 1. Geometry, design vector, per-animal scaling
%==========================================================================
nA = numel(mousedata);
if size(y, 1) ~= nA
    error('atlasVoxelCorrMap:ySize', 'y must have %d rows (one per animal).', nA);
end
if nA < 3
    error('atlasVoxelCorrMap:tooFewAnimals', 'Need at least 3 animals for a p-value.');
end
K = size(y, 2);

outsize = avsize;
if opts.Symmetrize
    if mod(avsize(3), 2) ~= 0
        error('atlasVoxelCorrMap:oddZ', 'Symmetrize requires an even avsize(3).');
    end
    outsize(3) = avsize(3)/2;
end
Ny = outsize(1); Nx = outsize(2); Nz = outsize(3);

if isempty(opts.Scale)
    scale = ones(nA, 1);
else
    scale = double(opts.Scale(:));
    assert(numel(scale) == nA, 'Scale must have one entry per animal.');
end

% Centre the design vector once: r = sum(xc.*yc) / (||xc|| * ||yc||)
yc  = y - mean(y, 1);
sy2 = sum(yc.^2, 1);                        % 1 x K
if any(sy2 == 0)
    error('atlasVoxelCorrMap:constantY', 'A column of y is constant.');
end

if isscalar(sigma), sigma = repmat(sigma, 1, 3); end
zrad = ceil(2*sigma(3));                    % imgaussfilt3 default FilterSize radius

%==========================================================================
% 2. Mask
%==========================================================================
haveMask = ~isempty(mask);
if haveMask
    assert(isequal(size(mask), outsize), ...
        'mask must be %dx%dx%d (the output volume size).', Ny, Nx, Nz);
    mask = logical(mask);
end

%==========================================================================
% 3. Pre-bin each animal's points by output z, so slab extraction is O(1)
%==========================================================================
pts = cell(nA, 1);
for a = 1:nA
    if isa(mousedata{a}, 'function_handle')
        pts{a} = mousedata{a};
    else
        pts{a} = prepPoints(mousedata{a}, outsize, opts.Symmetrize, scale(a));
    end
end

%==========================================================================
% 4. Slab thickness and GPU
%==========================================================================
if opts.SlabThickness > 0
    nzSlab = round(opts.SlabThickness);
else
    % live arrays per slab: Sx, Sxx, K x Sxy, one working volume (+ its buffer)
    bytesPerSlice = Ny*Nx*4;
    nArrays       = 3 + K;
    nzSlab        = floor(opts.MemoryBudgetGB*2^30/bytesPerSlice/nArrays) - 2*zrad;
end
nzSlab = max(1, min(Nz, nzSlab));

useGPU = opts.UseGPU;
if useGPU && ~hasGPU()
    warning('atlasVoxelCorrMap:noGPU', 'No GPU available, running on the CPU.');
    useGPU = false;
end

%==========================================================================
% 5. Slab loop
%==========================================================================
z0list = 1:nzSlab:Nz;
idxC   = cell(numel(z0list), 1);
rC     = cell(numel(z0list), 1);
muC    = cell(numel(z0list), 1);
sdC    = cell(numel(z0list), 1);
tstart = tic;

for s = 1:numel(z0list)
    z0    = z0list(s);
    z1    = min(Nz, z0 + nzSlab - 1);
    zl    = max(1,  z0 - zrad);             % buffered range actually built
    zh    = min(Nz, z1 + zrad);
    nzb   = zh - zl + 1;
    nzOut = z1 - z0 + 1;
    keepZ = (z0 - zl + 1) + (0:nzOut-1);    % crop the buffer back off

    % --- mask for this slab; skip empty slabs entirely --------------------
    if haveMask
        sel = mask(:,:,z0:z1);
        if ~any(sel(:))
            continue
        end
    else
        sel = true(Ny, Nx, nzOut);
    end

    % --- accumulators -----------------------------------------------------
    Sx  = zeros(Ny, Nx, nzOut, 'single');
    Sxx = zeros(Ny, Nx, nzOut, 'single');
    Sxy = zeros(Ny, Nx, nzOut, K, 'single');
    if useGPU
        Sx = gpuArray(Sx); Sxx = gpuArray(Sxx); Sxy = gpuArray(Sxy);
    end

    for a = 1:nA
        % build this animal's buffered slab ...
        if isa(pts{a}, 'function_handle')
            V = single(pts{a}(zl, zh)) * single(scale(a));
            assert(isequal(size(V), [Ny Nx nzb]), ...
                'Loader for animal %d returned %s, expected [%d %d %d].', ...
                a, mat2str(size(V)), Ny, Nx, nzb);
        else
            V = buildSlab(pts{a}, Ny, Nx, zl, zh);
        end
        if useGPU, V = gpuArray(V); end

        % ... smooth it (replicate padding only ever touches the buffer, or a
        % true volume edge where replicate is exactly what we want) ...
        V = imgaussfilt3(V, sigma, 'Padding', 'replicate');
        V = V(:,:,keepZ);

        % ... and accumulate the sufficient statistics.
        Sx  = Sx  + V;
        Sxx = Sxx + V.*V;
        for k = 1:K
            Sxy(:,:,:,k) = Sxy(:,:,:,k) + V * single(yc(a,k));
        end
    end

    % --- correlation over this slab ---------------------------------------
    ss   = max(Sxx - Sx.^2/nA, 0);          % sum of squared deviations
    dead = ss <= 0;                         % voxels with no across-animal variance
    if useGPU, selg = gpuArray(sel); else, selg = sel; end

    rSlab = zeros(nnz(sel), K, 'single');
    for k = 1:K
        rk       = Sxy(:,:,:,k) ./ sqrt(ss * single(sy2(k)));
        rk       = max(min(rk, 1), -1);
        rk(dead) = NaN;
        rSlab(:,k) = gather(rk(selg));
    end

    mu = Sx/nA;
    sd = sqrt(ss/(nA-1));

    idxC{s} = uint32(find(sel)) + uint32(z0-1)*uint32(Ny*Nx);
    rC{s}   = rSlab;
    muC{s}  = gather(mu(selg));
    sdC{s}  = gather(sd(selg));

    if opts.Verbose
        fprintf('atlasVoxelCorrMap: slab %d/%d (z %d-%d), %d voxels, %.1f s\n', ...
            s, numel(z0list), z0, z1, numel(idxC{s}), toc(tstart));
    end
    clear Sx Sxx Sxy ss dead selg V
end

%==========================================================================
% 6. Assemble, and turn r into p (cheap: only over the kept voxels)
%==========================================================================
out.idx     = cat(1, idxC{:});
out.r       = cat(1, rC{:});
% out.meanMap = cat(1, muC{:});
% out.stdMap  = cat(1, sdC{:});

df    = nA - 2;
r     = double(out.r);
t     = r .* sqrt(df ./ max(1 - r.^2, realmin));
pboth = betainc(df ./ (df + t.^2), df/2, 0.5);      % == 2*tcdf(-|t|, df)
switch opts.Tail
    case 'both'
        p = pboth;
    case 'right'
        p = pboth/2;  p(t <= 0) = 1 - pboth(t <= 0)/2;
    case 'left'
        p = pboth/2;  p(t >= 0) = 1 - pboth(t >= 0)/2;
end
p(~isfinite(r)) = NaN;
out.p = single(p);

out.n          = nA;
out.df         = df;
out.tail       = opts.Tail;
out.sigma      = sigma;
out.volumeSize = outsize;

if strcmp(opts.OutputFormat, 'full')
    out.rVol = nan(outsize, 'single');
    out.pVol = nan(outsize, 'single');
    out.rVol(out.idx) = out.r(:,1);
    out.pVol(out.idx) = out.p(:,1);
    if K > 1
        warning('atlasVoxelCorrMap:fullFirstColumn', ...
            'Full volumes were filled with the first column of y only.');
    end
end

if opts.Verbose
    fprintf('atlasVoxelCorrMap: %d voxels, %d animals, done in %.1f s\n', ...
        numel(out.idx), nA, toc(tstart));
end
end

%==========================================================================
function P = prepPoints(cellcoords, outsize, symflag, scaleval)
% Convert a point list into slice-sorted subscripts + weights, once per animal.

s = round(double(cellcoords(:, [2 1 3])));     % [row col page], as atlasCellMap
w = repmat(single(scaleval), size(s,1), 1);

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
P.rc  = uint32(s(ord,1) + (s(ord,2)-1)*outsize(1));   % linear index within a slice
P.z   = uint32(z);
P.w   = w(ord);
% points of slice k live at P.edge(k)+1 : P.edge(k+1)
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
