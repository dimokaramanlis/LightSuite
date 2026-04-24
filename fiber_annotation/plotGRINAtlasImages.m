function cf = plotGRINAtlasImages(slices_av, rvec_arr, radius_vox, depths_um, areaidx, namessub)
% PLOTGRINATLASIAMGES  Visualise atlas-region maps at multiple depths.
%
%   CF = plotGRINAtlasImages(SLICES_AV, RVEC_ARR, RADIUS_VOX, ...
%                            DEPTHS_UM, AREAIDX, NAMESSUB)
%
%   Each panel shows the colour-coded parcellation regions in a circular
%   cross-section at one depth below the lens bottom (0, 50, …, 250 µm).
%   A dashed circle marks the fiber boundary.  A shared colour-bar shows
%   the brain-area acronyms for all labelled regions.
%
%   Inputs
%   ------
%   SLICES_AV  – 1×Ndepths cell array of 2-D annotation arrays.
%   RVEC_ARR   – 1×Ndepths cell array of coordinate vectors (vox from centre).
%   RADIUS_VOX – fiber radius in atlas voxels (for the circle outline).
%   DEPTHS_UM  – 1×Ndepths vector of depth values in µm.
%   AREAIDX    – column vector of parcellation indices (from CSV).
%   NAMESSUB   – cell array of brain-area name strings, same order as AREAIDX.

    Ndepths = numel(depths_um);

    %------------------------------------------------------------------
    % Collect all unique parcellation IDs across all depth slices
    %------------------------------------------------------------------
    all_data = cat(3, slices_av{:});           % Ngrid × Ngrid × Ndepths

    [unareas, ~, iun_flat] = unique(all_data(:));
    iun_vol  = reshape(iun_flat, size(all_data));   % integer colour index

    store_areas  = unareas > 0;                % exclude outside-brain (ID 0)
    area_ids     = unareas(store_areas);       % parcellation IDs to label
    tick_idx     = find(store_areas);          % indices into unareas for caxis ticks

    Nareas = numel(area_ids);

    % Map each ID to an area name
    area_names = cell(Nareas, 1);
    for ii = 1:Nareas
        imatch = find(areaidx == area_ids(ii), 1);
        if ~isempty(imatch)
            area_names{ii} = namessub{imatch};
        else
            area_names{ii} = num2str(area_ids(ii));
        end
    end

    %------------------------------------------------------------------
    % Colourmap: one distinct colour per unique ID in iun_vol.
    % iun_vol is 1-indexed: value k → unareas(k).  unareas(1) is the
    % smallest ID, usually 0 (outside brain), which gets a near-black tone.
    %------------------------------------------------------------------
    Ncolors = numel(unareas);           % total unique values in iun_vol

    if Nareas > 0
        try
            area_cmap = cbrewer('qual', 'Set1', max(Nareas, 3));
            area_cmap = area_cmap(1:Nareas, :);
        catch
            area_cmap = lines(max(Nareas, 1));
        end
    else
        area_cmap = zeros(0, 3);
    end

    % Build full colourmap: one row per unique value in iun_vol.
    % Rows corresponding to IDs <= 0 get near-black; others get area colours.
    full_cmap = zeros(Ncolors, 3);
    area_row  = 0;
    for kk = 1:Ncolors
        if unareas(kk) <= 0
            full_cmap(kk, :) = [0.08 0.08 0.08];
        else
            area_row = area_row + 1;
            if area_row <= size(area_cmap, 1)
                full_cmap(kk, :) = area_cmap(area_row, :);
            else
                full_cmap(kk, :) = rand(1, 3);
            end
        end
    end

    %------------------------------------------------------------------
    % Figure layout
    %------------------------------------------------------------------
    panel_w = 220;
    fig_w   = min(1850, Ndepths * panel_w + 200);
    cf = figure('Name', 'GRIN Lens Atlas Regions', 'Color', 'w');
    cf.Position = [50, 80, fig_w, 420];

    pp = panel();
    pp.pack('h', Ndepths);
    pp.de.margin  = 4;
    pp.de.margintop = 16;
    pp.margin     = [8 8 8 32];
    pp.title('\bfGRIN Lens – Atlas Region Map');
    pp.fontname   = 'Arial';

    % Circle outline coordinates (in atlas voxels from centre)
    t_circ = linspace(0, 2*pi, 300);
    cx     = radius_vox * cos(t_circ);
    cy     = radius_vox * sin(t_circ);

    %------------------------------------------------------------------
    % Draw each depth panel
    %------------------------------------------------------------------
    clim_lo = 0.5;
    clim_hi = max(Ncolors, 1) + 0.5;

    for ii = 1:Ndepths
        rvec = rvec_arr{ii};
        iun  = iun_vol(:, :, ii);

        ax = pp(ii).select();
        imagesc(rvec, rvec, iun);
        axis(ax, 'equal', 'tight');
        ax.Visible = 'off';
        ax.Title.Visible = 'on';
        title(ax, sprintf('%d µm', depths_um(ii)), 'FontSize', 9, 'FontWeight', 'bold');

        ax.Colormap = full_cmap;
        clim(ax, [clim_lo, clim_hi]);

        % Fiber boundary circle
        line(ax, cx, cy, 'Color', 'w', 'LineWidth', 2, 'LineStyle', '--');
    end

    %------------------------------------------------------------------
    % Shared colourbar on the last panel
    %------------------------------------------------------------------
    if Nareas > 0
        last_ax = pp(Ndepths).select();
        cbar = colorbar(last_ax, 'eastoutside');
        cbar.Ticks      = tick_idx;
        cbar.TickLabels = area_names;
        cbar.FontSize   = 8;
        cbar.Label.String = 'Brain region';
    end
end
