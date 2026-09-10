function cf = plotCordCenterline(vol, mask, auto)
%PLOTCORDCENTERLINE Quality-control figure for the automatic cord centre line.
%
%   CF = PLOTCORDCENTERLINE(VOL, MASK, AUTO) draws the predictions returned by
%   AUTOCORDCENTERLINE on top of the data they came from, so a wrong automatic
%   answer is visible before anyone starts clicking in SPINAL_CORD_ALIGNER.
%
%   The figure has three parts:
%     * the segmented cross-sectional area along the cord, which is what the
%       rostrocaudal direction and the brain cut-off were decided from;
%     * the predicted centre (x and y) and angle along the cord;
%     * a row of sections spread over the cord, each with the predicted centre
%       and the predicted dorsoventral axis drawn on it. The green end of the
%       axis is the one the aligner treats as anterior.
%
%   The figure is created invisible and is meant to be printed to a file by the
%   caller.
%
%   See also AUTOCORDCENTERLINE, PREPARECORDSAMPLEFORREGISTRATION.

Nz    = size(vol, 3);
Nshow = 8;
ishow = round(linspace(0.02 * Nz, 0.98 * Nz, Nshow));
ishow = min(max(ishow, 1), Nz);

cf = figure('Visible', 'off');
cf.Position = [50 50 1500 750];

pp = panel();
pp.pack('v', {0.45 0.55});
pp(1).pack('h', 3);
pp(2).pack('h', Nshow);
pp.de.margin  = 12;
pp.margin     = [22 22 6 8];

slicevec = (1:Nz)';
%--------------------------------------------------------------------------
pp(1, 1).select();
plot(slicevec, auto.area, 'k-');
xlabel('Slice along the cord'); ylabel('Segmented area (px)');
title('Cross-section'); axis tight; grid on;
%--------------------------------------------------------------------------
pp(1, 2).select(); hold on;
plot(slicevec, auto.cen_x, 'b-', 'LineWidth', 1.2);
plot(slicevec, auto.cen_y, 'm-', 'LineWidth', 1.2);
plot(slicevec, auto.raw.cen(:, 1), 'b.', 'MarkerSize', 3);
plot(slicevec, auto.raw.cen(:, 2), 'm.', 'MarkerSize', 3);
xlabel('Slice along the cord'); ylabel('Position (px)');
legend({'x', 'y'}, 'Location', 'best');
title('Predicted centre'); axis tight; grid on;
%--------------------------------------------------------------------------
% The angle is drawn unwrapped along the cord, not wrapped to (-pi, pi]. A
% wrapped plot hides the one failure that matters here - net twist the images do
% not contain - behind jumps at the seam that look the same as real rotation.
thetaplot = getOr(auto, 'thetaunwr', unwrap(auto.theta));
pp(1, 3).select(); hold on;
plot(slicevec, auto.raw.theta, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
plot(slicevec, thetaplot, '-', 'Color', [0 0.6 0], 'LineWidth', 1.2);
xlabel('Slice along the cord'); ylabel('Angle (rad, unwrapped)');
title(sprintf('Predicted angle (net twist %2.0f deg)', ...
    rad2deg(thetaplot(end) - thetaplot(1))));
axis tight; grid on;
%--------------------------------------------------------------------------
rng(1);
isamp = randperm(numel(vol), min(numel(vol), 2e4));
vsamp = single(vol(isamp));
vmax  = quantile(vsamp, 0.999);
vmin  = quantile(vsamp, 0.001);

for ii = 1:Nshow
    islice = ishow(ii);
    pp(2, ii).select(); hold on;

    im = (single(vol(:, :, islice)) - vmin) / max(vmax - vmin, eps);
    imagesc(max(0, min(1, im)), [0 1]);
    contour(mask(:, :, islice), [0.5 0.5], 'Color', [0.3 0.6 1], 'LineWidth', 0.5);

    cx = auto.cen_x(islice);
    cy = auto.cen_y(islice);
    th = auto.theta(islice);
    r  = auto.rad(islice);
    if isnan(r) || r <= 0
        r = 0.2 * max(size(vol, [1 2]));
    end
    plot([cx - r*cos(th) cx + r*cos(th)], [cy - r*sin(th) cy + r*sin(th)], ...
        'y-', 'LineWidth', 1.2);
    plot(cx + r*cos(th), cy + r*sin(th), 'gs', 'MarkerSize', 5, 'LineWidth', 1);
    plot(cx, cy, 'c+', 'MarkerSize', 8, 'LineWidth', 1.2);

    ax = gca;
    ax.YDir = 'reverse'; ax.Colormap = gray; ax.Visible = 'off';
    ax.Title.Visible = 'on';
    axis equal; axis tight;
    title(sprintf('%d', islice));
end
%--------------------------------------------------------------------------
end
