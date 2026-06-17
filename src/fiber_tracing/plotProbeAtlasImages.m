function [cf, pp] = plotProbeAtlasImages(probe_ccf)
% PLOTPROBEATLASIMAGES  Visualise fitted Neuropixels probe trajectories.
%
%   [CF, PP] = plotProbeAtlasImages(PROBE_CCF)
%
%   Left panel  – 3-D Allen brain-grid wireframe with, for every probe, the
%                 annotated points and the fitted insertion→tip line, colour
%                 coded per probe.
%   Right panel – one vertical region column per probe showing the brain areas
%                 along the trajectory (colour-coded, labelled with acronyms)
%                 against depth from the brain surface in µm.
%
%   Input
%   -----
%   PROBE_CCF – struct array (one element per probe) produced by
%               annotateNeuropixelsProbes / fitProbeTrajectory.  Required
%               fields: points, trajectory_coords, trajectory_areas.  An
%               optional probe_number field is used for labelling.
%
%   Points and coordinates are in Allen CCF voxels with column order [AP DV ML]
%   and are mapped to the brain-grid axes as (AP, ML, DV) to match the atlas
%   orientation used elsewhere in LightSuite.

    Nprobes = numel(probe_ccf);
    if Nprobes == 0
        error('plotProbeAtlasImages: probe_ccf is empty.');
    end

    % Probe numbers (fall back to sequential order)
    if isfield(probe_ccf, 'probe_number')
        probe_nums = [probe_ccf.probe_number];
    else
        probe_nums = 1:Nprobes;
    end

    probe_cmap = neuropixels_probe_colors(max(probe_nums));

    %------------------------------------------------------------------
    % Brain-grid wireframe
    %------------------------------------------------------------------
    mf = which('brainGridData.npy');
    if isempty(mf)
        error('plotProbeAtlasImages: brainGridData.npy not found on the MATLAB path.');
    end
    bgdata = reduceBrainGrid(readNPY(mf), 6);

    %------------------------------------------------------------------
    % Figure layout
    %------------------------------------------------------------------
    fig_w = min(1600, 520 + Nprobes * 150);
    fig_h = 720;
    cf = figure('Name', 'Neuropixels Probe Trajectories', 'Color', 'w');
    cf.Position = [60, 80, fig_w, fig_h];

    pp = panel();
    pp.pack('h', {0.42 0.58});
    pp(2).pack('h', Nprobes);
    pp.fontname = 'Arial';
    pp.fontsize = 8;
    pp.de.margin = 10;
    pp.margin    = [14 14 6 6];
    pp(2).marginleft = 16;

    %------------------------------------------------------------------
    % Left: 3-D brain grid + probe points and fitted lines
    %------------------------------------------------------------------
    ax3d = pp(1).select();
    cla(ax3d);
    plotBrainGrid(bgdata, ax3d);
    hold(ax3d, 'on');

    for ii = 1:Nprobes
        col = probe_cmap(probe_nums(ii), :);

        pts = probe_ccf(ii).points;            % [AP DV ML]
        tc  = probe_ccf(ii).trajectory_coords; % 2×3 [AP DV ML]

        % Map [AP DV ML] -> brain-grid axes (AP, ML, DV)
        scatter3(ax3d, pts(:,1), pts(:,3), pts(:,2), 12, 'filled', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');
        plot3(ax3d, tc(:,1), tc(:,3), tc(:,2), '-', 'Color', col, 'LineWidth', 2.5);

        % Mark the tip (deepest point in DV)
        [~, itip] = max(tc(:,2));
        scatter3(ax3d, tc(itip,1), tc(itip,3), tc(itip,2), 45, col, 'filled', ...
            'MarkerEdgeColor', 'k');
    end
    title(ax3d, '3-D probe trajectories');

    %------------------------------------------------------------------
    % Right: per-probe region-depth columns
    %------------------------------------------------------------------
    max_depth = 0;
    for ii = 1:Nprobes
        max_depth = max(max_depth, max(probe_ccf(ii).trajectory_areas.trajectory_depth(:,2)));
    end
    if max_depth <= 0; max_depth = 1; end

    for ii = 1:Nprobes
        ax = pp(2, ii).select();
        hold(ax, 'on');

        ta = probe_ccf(ii).trajectory_areas;
        for r = 1:height(ta)
            d0 = ta.trajectory_depth(r, 1);
            d1 = ta.trajectory_depth(r, 2);
            c  = ta.color(r, :) / 255;
            patch(ax, 'XData', [0 1 1 0], 'YData', [d0 d0 d1 d1], ...
                'FaceColor', c, 'EdgeColor', 'none');

            % Label regions thick enough to be readable
            if (d1 - d0) > max_depth * 0.03
                txt_col = 'k';
                if mean(c) < 0.45; txt_col = 'w'; end
                text(ax, 0.5, (d0 + d1) / 2, ta.acronym{r}, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', 7, 'Color', txt_col, 'Interpreter', 'none');
            end
        end

        ax.YDir = 'reverse';
        ax.XLim = [0 1];
        ax.YLim = [0 max_depth * 1.02];
        ax.XTick = [];
        box(ax, 'on');
        col = probe_cmap(probe_nums(ii), :);
        title(ax, sprintf('Probe %d', probe_nums(ii)), 'Color', col);
        if ii == 1
            ylabel(ax, 'Depth from surface (µm)');
        else
            ax.YTickLabel = [];
        end
    end
end

%==========================================================================
function cmap = neuropixels_probe_colors(n)
% Distinct, stable colours for probes 1..n.  The first 9 match the group
% colours used by annotateNeuropixelsProbes so a probe keeps its colour
% between the annotation GUI and these plots.
    base = [ ...
        1.00 0.20 0.20;   % 1  red
        0.15 0.82 0.15;   % 2  green
        0.28 0.50 1.00;   % 3  blue
        1.00 0.55 0.00;   % 4  orange
        0.82 0.15 0.82;   % 5  magenta
        0.00 0.78 0.88;   % 6  cyan
        0.80 0.75 0.00;   % 7  yellow
        0.55 0.00 0.85;   % 8  purple
        0.00 0.65 0.45];  % 9  teal
    if n <= size(base, 1)
        cmap = base;
    else
        cmap = [base; lines(n - size(base, 1))];
    end
end
