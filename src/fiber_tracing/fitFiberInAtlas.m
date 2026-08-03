function [center, normal, info] = fitFiberInAtlas(pts, radius_vox, av, varargin)
% FITFIBERINATLAS  Fit a fiber / GRIN-lens axis to 3-D atlas-space points.
%
%   [CENTER, NORMAL]       = fitFiberInAtlas(PTS, RADIUS_VOX)
%   [CENTER, NORMAL, INFO] = fitFiberInAtlas(PTS, RADIUS_VOX, AV, 'Name', VALUE, ...)
%
%   PTS        – Nx3 clicked points in Allen CCF 10 µm voxel coords
%                [AP, DV, ML], collected across slices by annotateGRINLens.
%                They mark where the implant is – on its edge, inside it and
%                just around it.  They need not lie on any surface.
%   RADIUS_VOX – expected implant radius in 10 µm atlas voxels.  A guide: it
%                sets search scales and is reported next to the fitted value.
%   AV         – Allen CCF annotation volume, indexed AV(AP, DV, ML), 0
%                outside the brain.  Optional, but this is what supplies the
%                entry direction the axis is constrained against; without it
%                the constraint falls back to plain DV.
%
%   CENTER – 1x3 centre of the implant bottom face (tip) in atlas voxels.
%   NORMAL – 1x3 unit vector along the implant axis, pointing into the brain.
%   INFO   – diagnostics: .entry_normal .entry_point .tilt_deg .at_tilt_bound
%            .radius_fit .radius_ratio .inliers .axial_extent_vox
%            .depth_in_diameters .prior_weight ...
%
%   Name/value options
%   ------------------
%   'MaxTilt'     (45)      how far the axis may lean from the brain-surface
%                           normal at the entry, in degrees, per axis.
%   'Trim'        (0.95)    fraction of the clicks the fitted cylinder must
%                           cover; the rest are treated as strays.
%   'EntryPrior'  (0.50)    strength of the pull towards the entry normal on
%                           an annotation that carries no depth (see belowclo).
%                           0 fits the clicks alone.
%   'TipQuantile' (0.99)    quantile of the axial coordinate taken as the tip.
%   'OrientAxis'  ([0 1 0]) entry direction used when AV is not supplied.
%   'Verbose'     (true)    print a fit summary.
%
%   Why not PCA – and why a robust regression would not fix it either
%   -----------------------------------------------------------------
%   For a cloud filling a cylinder of length L and radius R the variance along
%   the axis is L^2/12 against R^2/4 in each radial direction, so the leading
%   eigenvector is the axis only when L > sqrt(3)*R, i.e. when the implant
%   enters deeper than ~0.87 diameters.  Below that PCA returns a *radial*
%   direction and the axis comes out ~90 deg wrong – the usual situation for a
%   short GRIN lens.
%
%   Swapping in a robust loss does not help, because the failure is not caused
%   by outliers: *any* fit that minimises the summed squared distance to the
%   line is PCA (that sum is trace(S) - n'*S*n, so minimising it maximises
%   n'*S*n), and a mean-square criterion always prefers the fat direction of a
%   short cloud.  What is needed is a criterion measuring how *wide* the cloud
%   is around the line rather than how far its points are on average.
%
%   Approach
%   --------
%   The axis is the line minimising the radius of the cylinder that covers a
%   fraction 'Trim' of the clicks.  That minimum sits at the true axis for any
%   aspect ratio – tilting a cylinder inside its covering cylinder can only
%   widen it – and it is robust by construction, since the outermost
%   1-'Trim' of the points (misclicks, spill around the implant) are ignored.
%
%   The width is measured as the mean of the ordered distances in the band
%   just below the 'Trim' quantile, not as the quantile itself (see
%   edgeWidth): everything above a plain quantile is discarded, so a tilted
%   axis can simply *shed* the end points that give the tilt away.  A useful
%   side effect is that the objective becomes an L-statistic, hence continuous
%   – the two clicks that swap rank are equidistant at the swap – so
%   finite-difference gradients behave.
%
%   That is solved as a constrained optimisation in four numbers (fmincon,
%   SQP): the two tilt angles of the direction away from the entry normal and
%   the two coordinates of the line's lateral position.  The ±'MaxTilt' limit
%   is a bound on the first two, enforced at every iterate rather than clipped
%   afterwards.  The objective is not convex in the tilt – a cloud can look
%   narrow around more than one direction – so the solve is started from the
%   best few points of a coarse scan of the tilt box.
%
%   The entry direction comes from the annotation: the brain surface is the
%   set of labelled voxels with an unlabelled neighbour, and the patch of it
%   nearest the clicks is locally flat, so its normal is the direction the
%   implant went in along.  Not every stretch of surface qualifies, though –
%   see the exposure test in brainEntryNormal.
%
%   For a long implant the constraint never binds and the clicks decide.  For
%   a short one the covering radius hardly changes with direction, so ±45° on
%   its own is not enough – inside that box the fit still wanders tens of
%   degrees on noise.  'EntryPrior' adds a quadratic pull back towards the
%   entry normal, weighted by how much the annotation can say about the lean
%   in the first place.  That weight is read off the *shape* of the click
%   cloud, rotation-invariantly: the ratio of its two largest scatter
%   eigenvalues is 1 for a ring or a blob, which fixes no direction at all,
%   and grows with the length of the annotated shaft.  A blob returns the
%   entry normal; a shaft annotated over a couple of diameters overrules the
%   prior completely.  In simulation that is worth 8.7° instead of 15° on a
%   real 25° lean, without letting a 0.4-diameter annotation run off to 50°.
%
%   Finally the point on the axis is slid along it to the 'TipQuantile' of the
%   clicks' axial coordinate, so CENTER lands on the bottom face of the
%   implant rather than in the middle of the annotation.  Only clicks the
%   fitted cylinder actually covers count, otherwise a single deep misclick
%   sets the tip.

    %------------------------------------------------------------------
    % 0. Options
    %------------------------------------------------------------------
    ip = inputParser;
    ip.FunctionName = 'fitFiberInAtlas';
    addRequired( ip, 'pts',        @(x) isnumeric(x) && ismatrix(x) && size(x,2) == 3);
    addRequired( ip, 'radius_vox', @(x) isnumeric(x) && isscalar(x) && x > 0);
    addOptional( ip, 'av',         [],       @(x) isempty(x) || isnumeric(x) || islogical(x));
    addParameter(ip, 'MaxTilt',     45,      @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 90);
    addParameter(ip, 'Trim',        0.95,    @(x) isnumeric(x) && isscalar(x) && x > 0.5 && x <= 1);
    addParameter(ip, 'EntryPrior',  0.50,    @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(ip, 'TipQuantile', 0.99,    @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(ip, 'OrientAxis',  [0 1 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(ip, 'Verbose',     true,    @(x) islogical(x) || isnumeric(x));
    if nargin < 3; av = []; end
    parse(ip, pts, radius_vox, av, varargin{:});
    opt = ip.Results;
    av  = opt.av;

    pts = double(pts);
    N   = size(pts, 1);
    if N < 4
        error('fitFiberInAtlas: need at least 4 points (got %d).', N);
    end
    R0 = double(radius_vox);

    %------------------------------------------------------------------
    % 1. Entry direction: inward normal of the brain surface where the
    %    implant goes in
    %------------------------------------------------------------------
    [n_entry, entry] = brainEntryNormal(av, pts, R0);
    if ~entry.found
        n_entry = opt.OrientAxis(:)' / norm(opt.OrientAxis);
    end

    %------------------------------------------------------------------
    % 2. Robust straight-line fit, constrained to the entry cone
    %------------------------------------------------------------------
    [normal, c_axis, R_cov, tilt, fitdiag] = constrainedAxis(pts, n_entry, ...
        opt.MaxTilt, opt.Trim, opt.EntryPrior);

    %------------------------------------------------------------------
    % 3. Slide the centre down the axis to the bottom face of the implant
    %------------------------------------------------------------------
    % The axis points into the brain, so the deepest clicks the fitted
    % cylinder actually covers mark the lens face.  Strays are kept out of
    % this, otherwise one misclick sets the tip.
    d   = perpDistance(pts, c_axis, normal);
    inl = d <= max(R_cov, 1);
    if nnz(inl) < 4; inl = true(N, 1); end

    t       = (pts - c_axis) * normal';
    tip_off = qvalue(t(inl), opt.TipQuantile);
    center  = c_axis + tip_off * normal;

    %------------------------------------------------------------------
    % 4. Diagnostics
    %------------------------------------------------------------------
    ang_entry = acosd(min(1, max(-1, normal * n_entry')));
    at_bound  = any(abs(abs(tilt) - opt.MaxTilt) < 0.25);
    extent    = max(t(inl)) - min(t(inl));

    info = struct();
    info.axis_centre        = c_axis;      % point on the fitted axis
    info.entry_normal       = n_entry;     % inward surface normal used
    info.entry_from_atlas   = entry.found; % false = fell back to OrientAxis
    info.entry_point        = entry.point; % surface voxel nearest the clicks
    info.entry_distance_vox = entry.dist;
    info.entry_flatness     = entry.flatness;
    info.entry_spread_deg   = entry.spread;  % how flat the surface is here
    info.tilt_deg           = tilt;        % [a b] lean from the entry normal
    info.tilt_total_deg     = ang_entry;
    info.at_tilt_bound      = at_bound;
    info.radius_fit         = R_cov;       % radius covering 'Trim' of the clicks
    info.radius_nominal     = R0;
    info.radius_ratio       = R_cov / R0;
    info.residuals          = d - R_cov;
    info.inliers            = inl;
    info.n_outliers         = nnz(~inl);
    info.axial_extent_vox   = extent;
    % How much depth the annotation covers, measured against the entry normal
    % rather than the fitted axis – the fitted axis tilts towards whatever
    % spread the clicks have, so measuring along it would let a flat
    % annotation report itself as a deep one.
    info.depth_in_diameters = fitdiag.depth_in_diameters;
    info.prior_weight       = fitdiag.prior_weight;   % 0 = clicks alone
    info.tip_offset_vox     = tip_off;
    info.n_points           = N;

    if opt.Verbose
        if entry.found
            fprintf('  fitFiberInAtlas: %d points | entry normal [%.2f %.2f %.2f] from %d surface voxels %.0f µm away (surface varies ±%.0f°)\n', ...
                N, n_entry, entry.n_vox, 10 * entry.dist, entry.spread);
        else
            fprintf('  fitFiberInAtlas: %d points | no brain surface found, entry normal defaults to [%.2f %.2f %.2f]\n', ...
                N, n_entry);
        end
        fprintf('    axis leans %.1f° from it (%.1f°, %.1f° per axis, limit ±%.0f°)\n', ...
            ang_entry, tilt(1), tilt(2), opt.MaxTilt);
        fprintf('    covering radius %.1f vox = %.0f µm diameter (%.2f x nominal) | %d stray point(s)\n', ...
            R_cov, 2 * 10 * R_cov, info.radius_ratio, info.n_outliers);
        fprintf('    annotated depth %.0f µm = %.1f diameters (surface prior weight %.2f)\n', ...
            10 * extent, info.depth_in_diameters, info.prior_weight);
    end

    if at_bound
        warning(['fitFiberInAtlas: the axis ran into the ±%.0f° limit, so the ' ...
                 'clicks want an implant far more oblique than the brain ' ...
                 'surface allows. Check that they all belong to one implant, ' ...
                 'or raise ''MaxTilt''.'], opt.MaxTilt);
    end
    if info.depth_in_diameters < 1
        warning(['fitFiberInAtlas: the clicks cover only %.1f implant diameters ' ...
                 'of depth. A flat annotation barely constrains the lean, so ' ...
                 'the axis mostly reflects the brain surface at the entry – ' ...
                 'click over a wider depth range, including the part of the ' ...
                 'implant above the surface if it is visible.'], ...
                 info.depth_in_diameters);
    end
end

%==========================================================================
% ROBUST CONSTRAINED LINE FIT
%==========================================================================

function [n, c, Rcov, tilt, diag] = constrainedAxis(pts, n0, maxtilt, q, prior)
% Line minimising the width of the click cloud around it, with the direction
% held within ±MAXTILT degrees of N0 in each of the two transverse axes.
%
% Posed as a constrained optimisation in four numbers – the two tilt angles
% and the two coordinates of the line's lateral position – so the ±MAXTILT
% limit is a plain bound that FMINCON enforces at every iterate, rather than
% something clipped afterwards.  The objective is an L-statistic (a weighted
% sum of order statistics), so although it has kinks where two clicks swap
% rank it never jumps: the two swapping points are equidistant at the swap.
% Finite-difference gradients cope with that provided the step is not
% vanishingly small, hence the explicit TypicalX.
    [u, v] = perpBasis(n0);
    c0     = geometricMedian(pts);      % robust reference point on the line

    %--- how much the clicks can say about the direction at all ----------
    % Deliberately rotation invariant: the ratio of the two largest scatter
    % eigenvalues is 1 for a ring or a disc of clicks, which carries no
    % direction information at all, and grows with the length of the
    % annotated shaft.  For a cloud filling a cylinder, lambda1/lambda2 is
    % 1 + (L/D)^2/0.75 once the shaft is long enough to dominate, which is
    % what the inversion below undoes.
    %
    % Measuring the depth along n0 instead – as an earlier version did – lets
    % a *wrong* entry normal foreshorten a deep annotation into an apparently
    % flat one and then over-weight its own prior on the strength of that,
    % which is exactly backwards.  It is what pulled fiber 1 of the YX047
    % dataset off: a 1.4 mm deep annotation reported itself as 0.61 diameters.
    % The strays have to come off first: a covariance is not robust, and three
    % misclicks at arm's length inflate the second eigenvalue enough to halve
    % the apparent elongation and so double the prior weight.
    dm     = sqrt(sum((pts - c0).^2, 2));
    core   = pts(dm <= qvalue(dm, q), :);
    if size(core, 1) < 4; core = pts; end
    lam    = sort(eig(cov(core)), 'descend');
    elong  = sqrt(lam(1) / max(lam(2), eps));
    aspect = sqrt(max(0, elong^2 - 1)) / 1.155;        % ~ length / diameter
    prior  = prior / (1 + aspect^2);

    diag = struct('prior_weight', prior, 'depth_in_diameters', aspect);

    %--- objective -------------------------------------------------------
    fun = @(x) edgeWidth(lineResiduals(x, pts, c0, n0, u, v), q) ...
               * (1 + prior * (x(1)^2 + x(2)^2) / maxtilt^2);

    %--- starting points -------------------------------------------------
    % The width is not convex in the tilt – a cloud can look narrow around
    % more than one direction – so the search starts from the best few points
    % of a coarse scan instead of a single guess.  Every candidate direction
    % is scored in one matrix product, so the scan is nearly free.
    g      = -maxtilt : maxtilt/4 : maxtilt;
    [A, B] = ndgrid(g, g);
    D = n0 + tand(A(:)) * u + tand(B(:)) * v;
    D = D ./ vecnorm(D, 2, 2);
    P = pts - c0;
    W = edgeWidth(sqrt(max(0, sum(P.^2, 2) - (P * D').^2)), q) ...
        .* (1 + prior * (A(:)'.^2 + B(:)'.^2) / maxtilt^2);
    [~, ord] = sort(W(:));
    nstart   = min(3, numel(ord));
    starts   = [A(ord(1:nstart)), B(ord(1:nstart)), zeros(nstart, 2)];

    %--- solve -----------------------------------------------------------
    w0 = max(qvalue(lineResiduals([0 0 0 0], pts, c0, n0, u, v), q), 1);
    lb = [-maxtilt, -maxtilt, -inf, -inf];
    ub = [ maxtilt,  maxtilt,  inf,  inf];

    xbest = starts(1, :);
    fbest = inf;
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', ...
            'Algorithm',                'sqp', ...
            'Display',                  'off', ...
            'FiniteDifferenceStepSize', 1e-2, ...
            'TypicalX',                 [maxtilt/4, maxtilt/4, w0, w0], ...
            'MaxFunctionEvaluations',   800, ...
            'StepTolerance',            1e-6);
        for s = 1:nstart
            [x, f] = fmincon(fun, starts(s,:), [], [], [], [], lb, ub, [], opts);
            if f < fbest; fbest = f; xbest = x; end
        end
    else
        % No Optimization Toolbox: Nelder-Mead on the same objective with the
        % tilt projected into the box.  Same answer whenever the optimum is
        % interior, which is the case whenever the fit is trustworthy anyway.
        box = @(x) [max(-maxtilt, min(maxtilt, x(1:2))), x(3:4)];
        for s = 1:nstart
            [x, f] = fminsearch(@(x) fun(box(x)), starts(s,:), ...
                optimset('Display', 'off', 'MaxFunEvals', 1200, 'TolX', 1e-4));
            x = box(x);
            if f < fbest; fbest = f; xbest = x; end
        end
    end

    %--- unpack ----------------------------------------------------------
    tilt = xbest(1:2);
    n    = n0 + tand(tilt(1)) * u + tand(tilt(2)) * v;
    n    = n / norm(n);
    c    = c0 + xbest(3) * u + xbest(4) * v;
    c    = c + mean((pts - c) * n') * n;    % park it at the middle of the cloud
    Rcov = qvalue(perpDistance(pts, c, n), q);
end

function d = lineResiduals(x, pts, c0, n0, u, v)
% Perpendicular distances from the clicks to the line described by X:
% X(1:2) tilt the direction away from N0 by that many degrees in the (n0,u)
% and (n0,v) planes – so a bound on them is literally "no more than 45° in
% either axis of entry" – and X(3:4) slide the line sideways.
    n = n0 + tand(x(1)) * u + tand(x(2)) * v;
    n = n / norm(n);
    P = pts - (c0 + x(3) * u + x(4) * v);
    d = sqrt(max(0, sum(P.^2, 2) - (P * n').^2));
end

function y = edgeWidth(D, q)
% How wide the cloud is around the line: the mean of the outermost distances
% in each column of D, up to the q-quantile.
%
% A single quantile is the obvious choice and is the wrong one.  Everything
% above it is discarded, so a tilted axis can simply *shed* the end points
% that give the tilt away - with 10% discarded the direction signal roughly
% halves.  The plain maximum keeps all of it but one misclick then sets the
% answer.  Averaging the band of order statistics just below the q-quantile
% keeps the outer envelope, which is where the direction information lives,
% while the 1-q above it absorbs the strays.
    n  = size(D, 1);
    k2 = min(n, max(1, round(q * n)));
    k1 = min(k2, max(1, round((q - 0.15) * n)));
    S  = sort(D, 1);
    y  = mean(S(k1:k2, :), 1);
end

function d = perpDistance(pts, c, n)
% Distance of every point from the line through C with direction N.
    P = pts - c;
    d = sqrt(max(0, sum(P.^2, 2) - (P * n').^2));
end

%==========================================================================
% BRAIN ENTRY DIRECTION, FROM THE ANNOTATION VOLUME
%==========================================================================

function [n0, entry] = brainEntryNormal(av, pts, radius_vox)
% Inward normal of the brain surface nearest the clicks.
%
% The annotation is 0 outside the brain and >0 inside, so the surface is
% simply the labelled voxels having an unlabelled 6-neighbour.  The patch of
% it around the entry is locally flat, and the normal of a plane *is* its
% smallest-variance direction – so an eigen-decomposition is the right tool
% here even though it is the wrong tool for the implant axis: a surface patch
% is a genuinely flat sheet, whereas the click cloud is not a genuine line.
%
% "Nearest" is not enough on its own, though – see the exposure test below.
    n0    = [0 1 0];
    entry = struct('found', false, 'point', [nan nan nan], 'dist', nan, ...
                   'n_vox', 0, 'flatness', nan, 'spread', nan);
    if isempty(av); return; end

    sz      = size(av);
    lo0     = min(pts, [], 1);
    hi0     = max(pts, [], 1);
    % Patch size is a real trade-off, settled against an independent estimate
    % (the gradient of a smoothed brain mask) at five cortical sites: 250 µm
    % is dominated by the voxel staircase and by whatever the surface happens
    % to be doing right at the entry - near the midline it reads the slope
    % into the interhemispheric fissure and comes back 40° off vertical - and
    % 1.2 mm averages across genuinely different parts of the surface (18°
    % errors at two of the five sites).  600 µm sat within 0.2-2.4° at all
    % five and is the default.
    patch_r = max(60, 2 * radius_vox);

    % A subsample of the clicks is plenty for locating the entry and keeps the
    % distance computation small for densely annotated implants.
    ps = pts(round(linspace(1, size(pts,1), min(size(pts,1), 100))), :);

    for pad = unique(min(250, max(30, 3 * radius_vox) * [1 2 4]))
        lo  = max([1 1 1], floor(lo0 - pad));
        hi  = min(sz,      ceil( hi0 + pad));
        sub = av(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) > 0;

        bnd = find(boundaryVoxels(sub));
        if isempty(bnd); continue; end      % box entirely inside the brain

        [i1, i2, i3] = ind2sub(size(sub), bnd);
        S = [i1, i2, i3] + lo - 1;          % global [AP DV ML] voxel coords

        % Distance of every surface voxel to the annotation, so the nearest
        % one can be taken as the entry.
        d2 = inf(size(S, 1), 1);
        for k = 1:size(ps, 1)
            d2 = min(d2, sum((S - ps(k,:)).^2, 2));
        end

        % Nearest surface first, but not every stretch of surface is an entry:
        % between the hemispheres, and around the ventricles, the annotation
        % has unlabelled sheets whose walls are "surface" too and often the
        % closest one to a medial implant.  Fitting a plane to those returns
        % the wall normal, ~90° from the real one.  A patch only qualifies if
        % there is open space beyond it, which is where the implant came from.
        avail = true(size(S, 1), 1);
        for attempt = 1:8
            dd = d2;  dd(~avail) = inf;
            [dmin, kmin] = min(dd);
            if ~isfinite(dmin); break; end

            entry_pt = S(kmin, :);
            near     = sum((S - entry_pt).^2, 2) <= patch_r^2;
            P        = S(near & avail, :);
            avail(near) = false;                 % this stretch is spent
            if size(P, 1) < 20; continue; end

            % Plane through the patch, refitted once without the 15% of
            % voxels sitting furthest off it.
            [n0, lam] = planeNormal(P);
            r  = abs((P - mean(P, 1)) * n0');
            kp = r <= qvalue(r, 0.85);
            if nnz(kp) >= 20
                P = P(kp, :);
                [n0, lam] = planeNormal(P);
            end

            % Orient it into the brain: march both ways from the patch voxels
            % themselves and keep whichever direction stays in labelled
            % tissue.  The march has to start *on* the surface – the patch
            % centroid sits below it wherever the surface curves, and from
            % there both directions look equally inside.
            Sm   = P(round(linspace(1, size(P,1), min(size(P,1), 200))), :);
            base = repmat(Sm, 10, 1);
            step = repelem((1:10)', size(Sm, 1), 1);
            walk = @(dir) sum(insideAnnotation(av, base + step * dir));
            if walk(n0) < walk(-n0); n0 = -n0; end

            % Open space beyond it?  Walk out from the patch centre: past a
            % real outer surface it is empty all the way, whereas a fissure
            % or ventricle wall runs back into tissue within a few voxels.
            m   = mean(P, 1);
            out = ~insideAnnotation(av, m - (6:36)' * n0);
            if mean(out) < 0.8; continue; end

            entry.found    = true;
            entry.point    = entry_pt;
            entry.dist     = sqrt(dmin);
            entry.n_vox    = size(P, 1);
            entry.flatness = lam(1) / max(lam(2), eps);   % 0 = perfectly flat
            entry.spread   = patchSpread(P, n0);
            return
        end
    end
    n0 = [0 1 0];       % nothing qualified – caller falls back to OrientAxis
end

function spread = patchSpread(S, n)
% How flat the surface actually is around the entry: split the patch into
% quadrants and see how far their normals disagree with the patch normal.
% A large value means the entry direction is only a rough guide – the surface
% is curving, or running down into a fissure – and the ±MaxTilt box is then
% doing more of the work than the normal itself.
    m = mean(S, 1);
    [u, v] = perpBasis(n);
    a = (S - m) * u' > 0;
    b = (S - m) * v' > 0;
    spread = 0;
    for qa = [false true]
        for qb = [false true]
            sel = a == qa & b == qb;
            if nnz(sel) >= 20
                nq = planeNormal(S(sel, :));
                spread = max(spread, acosd(min(1, abs(nq * n'))));
            end
        end
    end
end

function [n, lam] = planeNormal(S)
% Normal of the best-fitting plane through a set of voxels, plus the sorted
% eigenvalues so the caller can tell how plane-like the patch really was.
    C       = cov(S);
    [V, Dg] = eig((C + C') / 2);        % symmetric: eigenvalues ascending
    lam     = diag(Dg)';
    n       = V(:, 1)';
end

function bnd = boundaryVoxels(inside)
% Labelled voxels with at least one unlabelled 6-neighbour.
    bnd = false(size(inside));
    for dim = 1:3
        nd = size(inside, dim);
        if nd < 2; continue; end
        lo = repmat({':'}, 1, 3);  lo{dim} = 1:nd-1;
        hi = repmat({':'}, 1, 3);  hi{dim} = 2:nd;
        bnd(lo{:}) = bnd(lo{:}) | (inside(lo{:}) & ~inside(hi{:}));
        bnd(hi{:}) = bnd(hi{:}) | (inside(hi{:}) & ~inside(lo{:}));
    end
end

function ins = insideAnnotation(av, pts)
% True where a point falls inside labelled tissue; outside the volume counts
% as outside the brain.
    sz = size(av);
    ap = round(pts(:,1));
    dv = round(pts(:,2));
    ml = round(pts(:,3));
    ok = ap >= 1 & ap <= sz(1) & dv >= 1 & dv <= sz(2) & ml >= 1 & ml <= sz(3);
    ins = false(size(pts, 1), 1);
    if any(ok)
        ins(ok) = av(sub2ind(sz, ap(ok), dv(ok), ml(ok))) > 0;
    end
end

%==========================================================================
% SMALL HELPERS
%==========================================================================

function [u, v] = perpBasis(n)
% Orthonormal basis of the plane perpendicular to n.
    n   = n / norm(n);
    tmp = [1 0 0];
    if abs(n(1)) > 0.9; tmp = [0 1 0]; end
    u = cross(n, tmp);  u = u / norm(u);
    v = cross(n, u);    v = v / norm(v);
end

function c = geometricMedian(pts)
% Point minimising the summed distance to the cloud (Weiszfeld iteration) –
% a centre a handful of misclicks cannot drag away.
    c = median(pts, 1);
    for it = 1:100
        d  = max(sqrt(sum((pts - c).^2, 2)), 1e-6);
        cn = sum(pts ./ d, 1) / sum(1 ./ d);
        if norm(cn - c) < 1e-6; c = cn; return; end
        c = cn;
    end
end

function y = qvalue(x, q)
% q-quantile of a vector, MATLAB's convention, without the toolbox.
    y = colQuantile(x(:), q);
end

function y = colQuantile(X, q)
% q-quantile of every column of X, by linear interpolation between order
% statistics (the convention QUANTILE uses), vectorised over columns.
    n = size(X, 1);
    if n == 1; y = X; return; end
    S  = sort(X, 1);
    h  = min(max(q * n + 0.5, 1), n);
    lo = floor(h);
    hi = ceil(h);
    y  = S(lo, :) + (h - lo) * (S(hi, :) - S(lo, :));
end
