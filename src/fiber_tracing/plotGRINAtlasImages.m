function [cf, pp] = plotGRINAtlasImages(all_results,...
    areaidx, namessub, colorsub)
% PLOTGRINATLASIAMGES  Visualise atlas-region maps at multiple depths.
%
%   CF = plotGRINAtlasImages(SLICES_AV, RVEC_ARR, RADIUS_VOX, ...
%                            DEPTHS_UM, AREAIDX, NAMESSUB)
%   CF = plotGRINAtlasImages(..., FIBER_NAME)
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

    
    Ndepths = max(cellfun(@(x) numel(x.depths_um), all_results));
    Nfibers = numel(all_results);

    fiber_cmap = cbrewer('qual', 'Dark2', max(3, Nfibers));
    %------------------------------------------------------------------
    % Figure layout
    %------------------------------------------------------------------
    panel_w = 220;
    fig_w   = min(1850, Ndepths * panel_w + 200);
    fig_h   =  min(Nfibers*250, 1000);
    cf = figure('Name', ['Lens Atlas Regions'], 'Color', 'w');
    cf.Position = [50, 80, fig_w, fig_h];

    pp = panel();
    pp.pack('h', {0.4 0.6})
    pp(2).pack('v', Nfibers)
    for ifiber = 1:Nfibers 
        pp(2, ifiber).pack('h', {0.9 0.1});
        pp(2, ifiber,1).pack('h',Ndepths)
    end

    pp.fontname   = 'Arial';
    txtsize = 8;
    pp.fontsize = txtsize;
    pp.de.margin  = 1;
    pp(2).marginleft = 10;
    pp.margin     = [1 2 12 4];


    currax = pp(1).select(); cla;
    plotBrainGrid([], currax);
    hold on;
    
    for ii = 1:numel(all_results)
        ptsshow = all_results{ii}.atlas_pts;
        center = all_results{ii}.center_vox;
        w      = all_results{ii}.normal_vox;
        r      = all_results{ii}.radius_vox;
        
        % Ensure normal vector is a unit vector
        w = w(:)' / norm(w);
        
        % Build orthonormal basis for the plane perpendicular to the normal
        tmp = [1 0 0];
        if abs(dot(tmp, w)) > 0.9
            tmp = [0 1 0];
        end
        u = cross(w, tmp); u = u / norm(u);
        v = cross(w, u);   v = v / norm(v);
        
        % Determine cylinder span along the normal using projected points
        proj  = (ptsshow - center) * w';
        min_p = -(max(proj) + 100);
        max_p = 0;
        
        % Calculate 3D coordinates of the cylinder surface
        n_cyl = 30;
        theta = linspace(0, 2*pi, n_cyl);
        
        Cyl_AP = zeros(2, n_cyl);
        Cyl_DV = zeros(2, n_cyl);
        Cyl_ML = zeros(2, n_cyl);
        
        C1 = center + min_p * w;
        C2 = center + max_p * w;
        
        for k = 1:n_cyl
            offset = r * cos(theta(k)) * u + r * sin(theta(k)) * v;
            p1 = C1 + offset;
            p2 = C2 + offset;
            
            Cyl_AP(1, k) = p1(1); Cyl_DV(1, k) = p1(2); Cyl_ML(1, k) = p1(3);
            Cyl_AP(2, k) = p2(1); Cyl_DV(2, k) = p2(2); Cyl_ML(2, k) = p2(3);
        end
        
        % Plot the cylinder using DV(2), ML(3), AP(1) mapping
        surf(  Cyl_AP,Cyl_ML, Cyl_DV,'FaceColor', fiber_cmap(ii,:), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5);

        scatter3(ptsshow(:,2), ptsshow(:,3), ptsshow(:,1), 5, 'filled',...
            'MarkerFaceColor',fiber_cmap(ii,:),'MarkerEdgeColor','k')
    end

    for ifiber = 1:Nfibers
       %------------------------------------------------------------------
        % Collect all unique parcellation IDs across all depth slices
        %------------------------------------------------------------------
        slices_av = all_results{ifiber}.slices_av;
        radius_vox = all_results{ifiber}.radius_vox;
        rvec_arr   = all_results{ifiber}.rvec_arr;
        depths_um  = all_results{ifiber}.depths_um;
        all_data = cat(3, slices_av{:});           % Ngrid × Ngrid × Ndepths
    
        [unareas, ~, iun_flat] = unique(all_data(:));
        iun_vol  = reshape(iun_flat, size(all_data));   % integer colour index
    
        store_areas  = unareas > 0;                % exclude outside-brain (ID 0)
        area_ids     = unareas(store_areas);       % parcellation IDs to label
        tick_idx     = find(store_areas);          % indices into unareas for caxis ticks
    
        Nareas = numel(area_ids);
    
        % Map each ID to an area name
        area_names = cell(Nareas, 1);
        area_cmap  = nan(Nareas, 3);
        for ii = 1:Nareas
            imatch = find(areaidx == area_ids(ii), 1);
            if ~isempty(imatch)
                area_names{ii} = namessub{imatch};
                area_cmap(ii, :) = colorsub(imatch, :);
            else
                area_names{ii} = num2str(area_ids(ii));
                area_cmap(ii, :) = colorsub(num2str(area_ids(ii)), :);
            end
        end
        area_cmap = area_cmap/255;

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
    
       
        % Circle outline coordinates (in atlas voxels from centre)
        t_circ = linspace(0, 2*pi, 300);
        cx     = radius_vox * cos(t_circ);
        cy     = radius_vox * sin(t_circ);

        %------------------------------------------------------------------
        % Draw each depth panel
        %------------------------------------------------------------------
        clim_lo = 0.5;
        clim_hi = max(Ncolors, 1) + 0.5;
    
    
        axcolor = pp(2,ifiber,2).select(); 
        axcolor.Visible = 'off';
        positionbar = axcolor.Position;
        % positionbar(1) = positionbar(1) + positionbar(3)*0.8;
        positionbar(2) = positionbar(2)+positionbar(4)*0.25;
        positionbar(3) = positionbar(3)*0.2;
        positionbar(4) = positionbar(4)*0.5;
    % positionbar(4) = positionbar(4)*0.6;
    
        for ii = 1:Ndepths
            rvec = rvec_arr{ii};
            iun  = iun_vol(:, :, ii);
        
    
            ax = pp(2,ifiber, 1, ii).select();
            imagesc(rvec, rvec, iun);
    
            % 
            % gradim = imgradient(iun)~=0;
            % [row, col] = find(gradim);
            % line(rvec(col), rvec(row), 'Marker','.','Color','k','LineStyle','none')
    
            axis(ax, 'equal', 'tight');
            ax.Visible = 'off';
            ax.Title.Visible = 'on';
            title(ax, sprintf('%d µm', depths_um(ii)), 'FontSize', 9, 'FontWeight', 'bold');
    
            ax.Colormap = full_cmap;
            clim(ax, [clim_lo, clim_hi]);
    
            % Fiber boundary circle
            line(ax, cx, cy, 'Color', 'w', 'LineWidth', 2, 'LineStyle', '--');
            if ii == Ndepths
                cc = colorbar('Location', 'eastoutside');
                cc.Ticks      = tick_idx;
                cc.TickLabels = area_names;
                cc.FontSize   = txtsize-1;
                cc.Label.String = 'Brain region';
    
                cc.Box = 'off';
                cc.Position = positionbar;
    
            end
            if ii == 1
                ylabel(sprintf('Fiber %d', ifiber));
                yticks([])
                ax.YAxis.Visible = 'on';
                ax.YAxis.Color   = fiber_cmap(ifiber,:);
            end
        end
    end

  
end