function auto = autoCordCenterline(vol, mask, varargin)
%AUTOCORDCENTERLINE Automatic centre and dorsoventral axis of every cord section.
%
%   AUTO = AUTOCORDCENTERLINE(VOL, MASK) predicts, for every slice along the
%   long axis of the spinal cord volume VOL (dimension 3), where the cord sits
%   in the plane and how it is rotated. SPINAL_CORD_ALIGNER fits these
%   predictions as observations alongside the points the user clicks, dropping
%   them near any user point, so the user only has to correct the stretches
%   where the automatic answer is wrong.
%
%   Two quantities are estimated per slice:
%
%     centre  - the centroid of the segmented section. MASK is the binary
%               segmentation of the cord (Otsu threshold plus largest connected
%               component, see PREPARECORDSAMPLEFORREGISTRATION); the centroid
%               is the mean position of its pixels in that slice.
%
%     angle   - the direction of largest intensity variance inside the section.
%               The intensity-weighted covariance of the pixel positions is
%               formed with the within-section background subtracted, and its
%               leading eigenvector gives the axis.
%
%   An axis is only defined up to a half turn, and that has to be respected all
%   the way through or it invents twist that is not in the images: a section
%   lying near horizontal crosses the 0/pi seam for no physical reason, and
%   following it as an ordinary angle books that crossing as a 180 degree
%   rotation. So the axis is smoothed in doubled-angle space, where the seam does
%   not exist, and only then given a direction.
%
%   The direction comes from the third moment (skewness) of the intensity along
%   the axis - the angle points towards its heavier tail - but as a single vote
%   for the whole cord, weighted by how pronounced the asymmetry is on each
%   slice. Deciding it slice by slice does not work: the skewness is near zero
%   wherever the section is close to symmetric, so the direction would wander.
%
%   That leaves one global ambiguity, whether the predicted angle points
%   dorsally or ventrally, which depends on the stain and cannot be read off the
%   image at all. SPINAL_CORD_ALIGNER flips it for the whole cord with one key,
%   and any user-clicked slice overrides the prediction near it anyway.
%
%   AUTO = AUTOCORDCENTERLINE(..., 'lambda', L) sets the smoothing weight
%   (default 25). 'minpixels' (default 50) is the smallest section, in pixels,
%   that is still trusted to produce an estimate.
%
%   OUTPUT
%     auto - struct with fields
%              cen_x, cen_y - N x 1 predicted section centre (pixels, [x y]).
%              theta        - N x 1 predicted angle (rad), wrapped to (-pi, pi],
%                             in the convention of SPINAL_CORD_ALIGNER: the angle
%                             of the vector pointing from the posterior to the
%                             anterior mark.
%              thetaunwr    - the same angle, continuous along the cord. This is
%                             what anything fitting or interpolating the angle
%                             should use; unwrapping 'theta' again would only
%                             risk reintroducing the seam problem above.
%              rad          - N x 1 spread of the section along the axis, used
%                             only to draw the prediction.
%              area         - N x 1 number of segmented pixels per slice.
%              raw          - the unsmoothed per-slice estimates, for
%                             diagnostics: .axis is the bare axis angle, .theta
%                             the per-slice direction lifted onto the branch of
%                             the fit, .skew the asymmetry the vote used.
%              lambda       - the smoothing weight used.
%
%   See also PREPARECORDSAMPLEFORREGISTRATION, SPINAL_CORD_ALIGNER,
%   SMOOTHCORDSERIES.

%==========================================================================
p = inputParser;
addParameter(p, 'lambda',    25, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'minpixels', 50, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
lambda = p.Results.lambda;
minpix = p.Results.minpixels;
%==========================================================================
[Ny, Nx, Nz] = size(vol);
assert(isequal(size(mask), [Ny Nx Nz]), ...
    'autoCordCenterline:sizeMismatch', ...
    'The mask (%s) must have the same size as the volume (%s).', ...
    mat2str(size(mask)), mat2str([Ny Nx Nz]));

cen  = nan(Nz, 2);
phi  = nan(Nz, 1);   % major axis angle, direction not yet resolved
skw  = nan(Nz, 1);   % intensity skewness along the major axis
skw2 = nan(Nz, 1);   % and along the minor one
ani  = nan(Nz, 1);   % how well the two are told apart at all
rad  = nan(Nz, 1);
area = zeros(Nz, 1);

[colgrid, rowgrid] = meshgrid(single(1:Nx), single(1:Ny));
%==========================================================================
for islice = 1:Nz
    m         = mask(:, :, islice);
    area(islice) = nnz(m);
    if area(islice) < minpix
        continue
    end
    idx = find(m);
    xs  = colgrid(idx);
    ys  = rowgrid(idx);

    % --- centre of the section: centroid of the segmentation ---
    cen(islice, :) = [mean(xs) mean(ys)];

    % --- direction of intensity variance inside the section ---
    scurr = single(medfilt2(vol(:, :, islice), [3 3], 'symmetric'));
    w     = scurr(idx);
    w     = w - quantile(w, 0.05);   % background level *inside* the cord
    w(w < 0) = 0;
    sw    = sum(w);
    if sw <= 0
        continue
    end

    mx = sum(w .* xs) / sw;
    my = sum(w .* ys) / sw;
    dx = xs - mx;
    dy = ys - my;

    cxx = sum(w .* dx .* dx) / sw;
    cyy = sum(w .* dy .* dy) / sw;
    cxy = sum(w .* dx .* dy) / sw;

    [V, D]     = eig(double([cxx cxy; cxy cyy]));
    evals      = diag(D);
    [ev, imax] = max(evals);
    imin       = 3 - imax;
    u          = V(:, imax);

    % How much the two eigenvalues differ is how much of an axis there is at
    % all: on a section whose intensity spreads equally in every direction the
    % principal direction is arbitrary and can point anywhere without the image
    % changing. Recording that lets the fit ignore those sections later.
    ani(islice) = (ev - evals(imin)) / max(ev + evals(imin), eps);

    % skewness of the intensity fixes which way an axis points, and it is needed
    % for both principal directions, since which of them the cord is actually
    % followed along is only decided later
    t1 = dx * V(1, imax) + dy * V(2, imax);
    t2 = dx * V(1, imin) + dy * V(2, imin);
    skw(islice)  = weightedSkew(w, t1, sw);
    skw2(islice) = weightedSkew(w, t2, sw);

    phi(islice) = atan2(u(2), u(1));
    rad(islice) = sqrt(max(ev, 0));
end
%==========================================================================
% smooth the centre
icen = find(~isnan(cen(:, 1)));
if numel(icen) < 2
    error('autoCordCenterline:noSections', ...
        ['Only %d slices held a segmented section of at least %d pixels, which ' ...
         'is not enough to predict a centre line. Check the segmentation of the ' ...
         'registration channel.'], numel(icen), minpix);
end

auto       = struct();
auto.cen_x = smoothCordSeries(Nz, icen, cen(icen, 1), lambda);
auto.cen_y = smoothCordSeries(Nz, icen, cen(icen, 2), lambda);
%==========================================================================
% smooth the axis, in doubled-angle space
%
% phi and phi+pi describe the same axis, so the 0/pi seam is not a rotation and
% the cord crosses it wherever the section happens to lie near horizontal.
% Smoothing phi as an ordinary angle turns every one of those crossings into a
% half turn, and along a cord of a couple of thousand slices they add up into
% hundreds of degrees of twist that is not in the images at all. Doubling the
% angle removes the ambiguity - 2*phi is a genuine circular quantity for an axis
% - and smoothing its cosine and sine as vectors cannot drift the way
% unwrapping can.
%
% The observations are also weighted by how well each section defines an axis
% at all. The angular precision of a principal direction goes with the squared
% relative gap between the eigenvalues, so that is the weight: a section whose
% intensity spreads equally in every direction contributes almost nothing, and
% the fit coasts through it on the neighbours it can actually see. Without this
% a stretch of near-symmetric sections - the enlargement where a cord runs into
% the brainstem, typically - makes the prediction wander through a half turn or
% more that is nowhere in the images.
iax = find(~isnan(phi));
if numel(iax) > 1
    % Which of the two principal directions to follow is decided by continuity
    % along the cord, not by taking the larger one slice by slice. A cord
    % section is taller than it is wide in some places and wider than tall in
    % others, so "the direction of largest variance" flips by 90 degrees
    % wherever the aspect ratio crosses over - twice in a typical sample, once
    % into an enlargement and once out of it - and the fit then follows a
    % quarter turn that is not a twist of the cord at all. Tracking one physical
    % axis through those crossings keeps it pointing the same way.
    strack   = resolveAxisTrack(phi, ani, Nz);
    phitrack = phi + (pi/2) * (strack < 0);
    skwtrack = skw;
    skwtrack(strack < 0) = skw2(strack < 0);

    psi   = 2 * phitrack(iax);
    wax   = ani(iax).^2;
    % normalised to average 1, so LAMBDA keeps meaning the same thing as it does
    % for the unweighted series and only the relative weighting bites
    wax   = wax / max(mean(wax, 'omitnan'), eps);
    cpsi  = smoothCordSeries(Nz, iax, cos(psi), lambda, wax);
    spsi  = smoothCordSeries(Nz, iax, sin(psi), lambda, wax);

    % back to an axis, made continuous along the cord. Unwrapping here is safe
    % where unwrapping the raw estimates was not: the smoothed doubled angle
    % moves by a tiny amount from slice to slice by construction.
    phism = unwrap(atan2(spsi, cpsi)) / 2;

    %----------------------------------------------------------------------
    % Which end of the axis is anterior is decided once, for the whole cord.
    % The intensity skewness says so only weakly - it is near zero wherever the
    % section is close to symmetric - so a per-slice decision is mostly noise
    % and lets the direction wander. A vote weighted by how pronounced the
    % asymmetry actually is on each slice is far steadier, and a wrong global
    % choice is one key press in SPINAL_CORD_ALIGNER.
    dirpref = phitrack(iax) + pi * (skwtrack(iax) < 0);
    agrees  = sign(cos(dirpref - phism(iax)));
    votewt  = abs(skwtrack(iax)) .* wax;   % also discount sections with no axis
    score   = sum(votewt .* agrees, 'omitnan');
    if score < 0
        phism = phism + pi;
    end
    fprintf('Direction of the axis chosen by %2.0f%% of the weighted asymmetry.\n', ...
        100 * (0.5 + 0.5 * abs(score) / max(sum(votewt, 'omitnan'), eps)));

    auto.thetaunwr = phism;                    % continuous along the cord
    auto.theta     = wrapToPiLocal(phism);
    % the per-slice estimate, lifted onto the same branch as the fit, so raw and
    % fitted can be compared without the seam getting in the way
    rawlift        = nan(Nz, 1);
    rawlift(iax)   = phism(iax) + wrapToPiLocal(dirpref - phism(iax));
else
    auto.thetaunwr = zeros(Nz, 1);
    auto.theta     = zeros(Nz, 1);
    rawlift        = nan(Nz, 1);
    warning('autoCordCenterline:noAngle', ...
        'No slice produced an intensity axis; the predicted angle is 0 everywhere.');
end
%==========================================================================
irad = find(~isnan(rad));
if numel(irad) > 1
    auto.rad = smoothCordSeries(Nz, irad, rad(irad), lambda);
else
    auto.rad = repmat(mean(rad, 'omitnan'), Nz, 1);
end

auto.area   = area;
auto.lambda = lambda;
auto.raw    = struct('cen', cen, 'theta', rawlift, 'axis', phi, ...
                     'rad', rad, 'skew', skw, 'aniso', ani);
%==========================================================================
fprintf('Predicted the cord centre line from %d/%d sections. Net twist %2.0f deg.\n', ...
    numel(icen), Nz, rad2deg(auto.thetaunwr(end) - auto.thetaunwr(1)));
%==========================================================================
end

% =========================================================================
function sk = weightedSkew(w, t, sw)
%WEIGHTEDSKEW Intensity skewness along one direction, 0 when there is no spread.
sg = sqrt(sum(w .* t.^2) / sw);
if sg > 0
    sk = sum(w .* t.^3) / sw / sg^3;
else
    sk = 0;
end
end

% -------------------------------------------------------------------------
function strack = resolveAxisTrack(phi, ani, Nz)
%RESOLVEAXISTRACK Follow one principal direction along the cord, not the larger.
%   Returns +1 where the major axis continues the track and -1 where the minor
%   one does. The two candidates are 90 degrees apart, which is a half turn in
%   doubled-angle space, so this is a sign choice on the doubled-angle unit
%   vector of each slice.
%
%   The track is grown outwards from the most clearly oriented slice, carrying a
%   running reference that each new slice is compared against. The reference is
%   updated in proportion to how well that slice defines an axis, so a stretch
%   of near-circular sections neither flips the track nor destroys it: the
%   reference simply coasts through unchanged and comes out the far side still
%   pointing the way it went in.

alpha  = 0.15;
strack = ones(Nz, 1);

uvec  = [cos(2*phi), sin(2*phi)];
igood = find(all(isfinite(uvec), 2));
if numel(igood) < 2
    return
end

wn = ani;
wn(~isfinite(wn)) = 0;
wn = wn.^2 / max(median(wn(igood).^2), eps);
wn = min(wn, 4);              % one very anisotropic slice must not take over

[~, ianchor] = max(wn);
if ~any(igood == ianchor)
    ianchor = igood(1);
end

for pass = 1:2
    if pass == 1
        seq = ianchor:Nz;
    else
        seq = ianchor:-1:1;
    end
    ref = uvec(ianchor, :);
    for k = seq
        if ~all(isfinite(uvec(k, :)))
            continue
        end
        if dot(ref, uvec(k, :)) < 0
            strack(k) = -1;
        else
            strack(k) = 1;
        end
        beta = min(alpha * wn(k), 1);
        ref  = (1 - beta) * ref + beta * strack(k) * uvec(k, :);
        nref = norm(ref);
        if nref > eps
            ref = ref / nref;
        end
    end
end

end
