function [cf, pp] = plotGRINIntensityImages(all_results, intensity_results, channames)
% PLOTGRININTENSITYIMAGES  Visualise fluorescence intensity profiles at
%   multiple depths through each GRIN fiber in atlas space.
%
%   Layout
%   ------
%   Left  (35%) – 3-D brain-grid plot with colour-coded fiber cylinders.
%   Right (65%) – one row-group per fiber:
%       Panels (82%) – one row per channel, one column per depth.
%                      Each panel shows the grayscale intensity cross-section
%                      with atlas region boundaries overlaid in white.
%                      The within-circle median intensity is printed in the
%                      bottom-left corner of each panel.
%       Profile (18%) – line plot of median intensity vs depth for every
%                       channel of that fiber.
%
%   Contrast limits are shared across all fibers for the same channel
%   (1st–99th percentile of all nonzero in-circle voxels pooled across
%   every fiber and every depth).

    Nfibers   = numel(all_results);
    Nchannels = numel(channames);
    Ndepths   = numel(intensity_results{1}.depths_um);

    fiber_cmap = cbrewer('qual', 'Dark2', max(3, Nfibers));

    % Per-channel colormaps for intensity images
    chan_cmaps  = {'gray', 'hot', 'winter', 'copper', 'bone'};
    % Per-channel line colours for depth-profile plots (distinct, colormap-agnostic)
    chan_colors = [0.45 0.45 0.45;   % gray
                   0.85 0.25 0.05;   % red-orange
                   0.05 0.45 0.80;   % blue
                   0.10 0.65 0.20;   % green
                   0.60 0.05 0.80];  % purple

    %----------------------------------------------------------------------
    % Per-channel contrast limits – pooled across ALL fibers and depths
    %----------------------------------------------------------------------
    chan_clim = zeros(Nchannels, 2);
    for ichan = 1:Nchannels
        all_vals = [];
        for ifiber = 1:Nfibers
            for id = 1:Ndepths
                v = intensity_results{ifiber}.slices_int{id}(:,:,ichan);
                all_vals = [all_vals; v(v > 0)]; %#ok<AGROW>
            end
        end
        if isempty(all_vals)
            chan_clim(ichan,:) = [0, 1];
        else
            lo = double(quantile(all_vals, 0.01));
            hi = double(quantile(all_vals, 0.99));
            if hi <= lo; hi = lo + 1; end
            chan_clim(ichan,:) = [lo, hi];
        end
    end

    %----------------------------------------------------------------------
    % Figure layout
    %----------------------------------------------------------------------
    panel_w = 220;
    fig_w   = min(1850, Ndepths * panel_w + 350);
    fig_h   = min(Nfibers * Nchannels * 150 + 80, 1200);
    cf = figure('Name', 'Lens Intensity Profiles', 'Color', 'w');
    cf.Position = [100, 60, fig_w, fig_h];

    %%
    pp = panel();
    pp.pack('h', {0.35 0.65});
    pp(2).pack('v', Nfibers);
    for ifiber = 1:Nfibers
        pp(2, ifiber).pack('h', {0.82 0.18});
        pp(2, ifiber, 1).pack('v', Nchannels);
        for ichan = 1:Nchannels
            pp(2, ifiber, 1, ichan).pack('h', Ndepths);
        end
    end

    pp.fontname  = 'Arial';
    pp.fontsize  = 8;
    pp.de.margin = 1;
    pp(2).marginleft = 10;
    for ifiber = 1:Nfibers
        pp(2, ifiber,2).marginleft = 15;
        pp(2,ifiber).margintop = 15;
    end
    pp.margin    = [1 10 3 4];

    %----------------------------------------------------------------------
    % Left panel: 3-D brain grid + fiber cylinders
    %----------------------------------------------------------------------
    currax = pp(1).select();
    cla(currax);
    plotBrainGrid([], currax);
    hold(currax, 'on');

    for ii = 1:Nfibers
        ptsshow = all_results{ii}.atlas_pts;
        center  = all_results{ii}.center_vox;
        w       = all_results{ii}.normal_vox;
        r       = all_results{ii}.radius_vox;

        w = w(:)' / norm(w);
        tmp = [1 0 0];
        if abs(dot(tmp, w)) > 0.9; tmp = [0 1 0]; end
        u  = cross(w, tmp); u  = u  / norm(u);
        vv = cross(w, u);   vv = vv / norm(vv);

        proj  = (ptsshow - center) * w';
        min_p = -(max(proj) + 100);

        n_cyl  = 30;
        theta  = linspace(0, 2*pi, n_cyl);
        Cyl_AP = zeros(2, n_cyl);
        Cyl_DV = zeros(2, n_cyl);
        Cyl_ML = zeros(2, n_cyl);
        C1 = center + min_p * w;
        C2 = center;
        for k = 1:n_cyl
            off = r * cos(theta(k)) * u + r * sin(theta(k)) * vv;
            p1  = C1 + off;  p2 = C2 + off;
            Cyl_AP(1,k) = p1(1); Cyl_DV(1,k) = p1(2); Cyl_ML(1,k) = p1(3);
            Cyl_AP(2,k) = p2(1); Cyl_DV(2,k) = p2(2); Cyl_ML(2,k) = p2(3);
        end
        surf(currax, Cyl_AP, Cyl_ML, Cyl_DV, ...
            'FaceColor', fiber_cmap(ii,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        scatter3(currax, ptsshow(:,2), ptsshow(:,3), ptsshow(:,1), 5, 'filled', ...
            'MarkerFaceColor', fiber_cmap(ii,:), 'MarkerEdgeColor', 'k');
    end

    %----------------------------------------------------------------------
    % Right panels: intensity images and depth profiles
    %----------------------------------------------------------------------
    t_circ = linspace(0, 2*pi, 300);

    minint = inf;
    maxint = 0;
    for ifiber = 1:Nfibers
        int_res    = intensity_results{ifiber};
        currmin = min(int_res.median_intensity,[],'all');
        currmax = max(int_res.median_intensity,[],'all');
        maxint = max(maxint,currmax);
        minint = min(minint,currmin);
    end
    minint = floor(minint/50)*50;
    maxint = ceil(maxint/50)*50;

    for ifiber = 1:Nfibers
        int_res    = intensity_results{ifiber};
        depths_um  = int_res.depths_um;
        radius_vox = int_res.radius_vox;
        med_int    = int_res.median_intensity;   % [Ndepths × Nchannels]
        cx = radius_vox * cos(t_circ);
        cy = radius_vox * sin(t_circ);

        % ---- intensity panels ----
        for ichan = 1:Nchannels
            clim_lo   = chan_clim(ichan, 1);
            clim_hi   = chan_clim(ichan, 2);
            cmap_name = chan_cmaps{mod(ichan - 1, numel(chan_cmaps)) + 1};

            for id = 1:Ndepths
                rvec       = int_res.rvec_arr{id};
                slice_data = double(int_res.slices_int{id}(:,:,ichan));

                ax = pp(2, ifiber, 1, ichan, id).select();
                imagesc(ax, rvec, rvec, slice_data, [clim_lo, clim_hi]);
                colormap(ax, cmap_name);
                axis(ax, 'equal', 'tight');
                hold(ax, 'on');
                ax.Visible = 'off';

                % Depth label on the very first row only
                if ifiber == 1 && ichan == 1
                    ax.Title.Visible = 'on';
                    title(ax, sprintf('%d µm', depths_um(id)), ...
                        'FontSize', 9, 'FontWeight', 'bold');
                end

                % Atlas region boundary overlay
                ann = double(all_results{ifiber}.slices_av{id});
                bnd = false(size(ann));
                bnd(1:end-1, :) = bnd(1:end-1,:) | (ann(1:end-1,:) ~= ann(2:end,:));
                bnd(2:end,   :) = bnd(2:end,:)   | (ann(1:end-1,:) ~= ann(2:end,:));
                bnd(:, 1:end-1) = bnd(:,1:end-1) | (ann(:,1:end-1) ~= ann(:,2:end));
                bnd(:, 2:end  ) = bnd(:,2:end)   | (ann(:,1:end-1) ~= ann(:,2:end));
                bnd = bnd & (ann > 0);
                h_bnd = imagesc(ax, rvec, rvec, ones([size(bnd), 3]));
                h_bnd.AlphaData = double(bnd) * 0.55;

                % Fiber boundary circle
                line(ax, cx, cy, 'Color', fiber_cmap(ifiber,:), ...
                    'LineWidth', 1.5, 'LineStyle', '--');

                % Median intensity annotation (bottom-left corner)
                text(ax, 0.05, 0.05, sprintf('%.0f', med_int(id, ichan)), ...
                    'Units', 'normalized', 'FontSize', 7, 'Color', 'w', ...
                    'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

                % Row label on first depth column
                if id == 1
                    ylabel(ax, sprintf('F%d / %s', ifiber, channames{ichan}), ...
                        'FontSize', 8);
                    yticks(ax, []);
                    ax.YAxis.Visible = 'on';
                    ax.YAxis.Color   = fiber_cmap(ifiber,:);
                end
            end
        end

        % ---- depth profile (right column) ----
        ax_p = pp(2, ifiber, 2).select();
        hold(ax_p, 'on');
        for ichan = 1:Nchannels
            col = chan_colors(mod(ichan - 1, size(chan_colors,1)) + 1, :);
            plot(ax_p, med_int(:, ichan), depths_um, ...
                '-o', 'Color', col, 'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', channames{ichan});
        end
        ax_p.YDir    = 'reverse';
        ax_p.YLim    = [depths_um(1) - 10, depths_um(end) + 10];
        ax_p.XLim    =[minint maxint];
        ax_p.FontSize = 7;
        xlabel(ax_p, 'Median', 'FontSize', 7);
        ylabel(ax_p, 'Depth (µm)', 'FontSize', 7);
        title(ax_p, sprintf('Fiber %d', ifiber), 'FontSize', 8);
        box(ax_p, 'off');
        if Nchannels > 1
            legend(ax_p, 'Location', 'best', 'FontSize', 6, 'Box', 'off');
        end
    end
    %%
end
