function spinal_cord_aligner(opts)
% SPINAL_ALIGNER_OPTIMIZED 
%
%   CONTROLS:
%       Left Click   : Set ANTERIOR/FRONT (Green)
%       Right Click  : Set POSTERIOR/BACK (Red)
%       Middle Click : Set CENTER CANAL (Blue)
%       'p'          : Toggle Prediction visibility
%       'c'          : Clear all points on current slice
%       'x'          : Delete point nearest to cursor
%       Space        : Jump to largest gap
%       's'          : Save
%       [New] Click on side plots: Jump to specific slice

    % -- 1. Setup & Input Check --
    if nargin < 1 || ~isfield(opts, 'regvol')
        error('Input must be a struct with .regvol and .lsfolder fields.');
    end
    
    gui_data = struct();
    gui_data.save_path = opts.lsfolder;
    if ~isfolder(gui_data.save_path)
        mkdir(gui_data.save_path);
    end

    % Parse Volume
    [gui_data.Rows, gui_data.Cols, gui_data.Nslices] = size(opts.regvol);

    % Contrast Normalization
    rng(1);
    samp_idx = randperm(numel(opts.regvol), min(numel(opts.regvol), 2e4));
    v_samp   = single(opts.regvol(samp_idx));
    vmax     = quantile(v_samp, 0.999);
    vmin     = quantile(v_samp, 0.001);
    v_float  = (single(opts.regvol) - vmin) / (vmax - vmin);
    v_float  = max(0, min(1, v_float));
    gui_data.disp_vol = uint8(255 * v_float);
   
   % -- Data Structures --
    gui_data.user.cen = nan(gui_data.Nslices, 2); 
    gui_data.user.ant = nan(gui_data.Nslices, 2); 
    gui_data.user.pos = nan(gui_data.Nslices, 2);
    
    % Derived Observations (for plotting raw user points)
    gui_data.obs.x = nan(gui_data.Nslices, 1);
    gui_data.obs.y = nan(gui_data.Nslices, 1);
    gui_data.obs.th = nan(gui_data.Nslices, 1);
    
    % Solver outputs
    gui_data.fit.x     = nan(gui_data.Nslices, 1);
    gui_data.fit.y     = nan(gui_data.Nslices, 1);
    gui_data.fit.theta = nan(gui_data.Nslices, 1);
    gui_data.fit.rad   = nan(gui_data.Nslices, 1);
    
    % Settings
    gui_data.curr_slice = 1;
    gui_data.show_pred  = true; 
    
    % Default Regularization (overwritten if loaded)
    gui_data.lambda_pos = 5000;
    gui_data.lambda_ang = 5000;
    
    % -- 2. Check for Previous Save --
    save_file = fullfile(gui_data.save_path, 'spinal_alignment_opt.mat');
    if exist(save_file, 'file')
        answer = questdlg('Previous alignment found. Load it?', 'Load Data', 'Yes', 'No', 'Yes');
        if strcmp(answer, 'Yes')
            try
                loaded = load(save_file);
                d = loaded.align_out;
                if size(d.user_cen, 1) == gui_data.Nslices
                    gui_data.user.cen = d.user_cen;
                    gui_data.user.ant = d.user_ant;
                    gui_data.user.pos = d.user_pos;
                    
                    if isfield(d, 'lambda_pos')
                        gui_data.lambda_pos = d.lambda_pos;
                    end
                    if isfield(d, 'lambda_ang')
                        gui_data.lambda_ang = d.lambda_ang;
                    end
                    
                    fprintf('Data loaded (Reg: Pos=%.0f, Ang=%.0f).\n', gui_data.lambda_pos, gui_data.lambda_ang);
                    gui_data = run_optimizer(gui_data);
                end
            catch
                warning('Load failed.');
            end
        end
    end
    
    % -- 3. Initialize GUI --
    screen = get(0, 'ScreenSize');
    % Make the figure wider to accommodate side-by-side panels
    fig_dim_w = min(1200, screen(3)*0.95); 
    fig_dim_h = min(800, screen(4)*0.85);
    
    gui_fig = figure('Name', 'Spinal Cord Optimizer', ...
        'NumberTitle', 'off', ...
        'Position', [(screen(3)-fig_dim_w)/2, (screen(4)-fig_dim_h)/2, fig_dim_w, fig_dim_h], ...
        'Color', 'w', ...
        'MenuBar', 'none', ...
        'WindowScrollWheelFcn', @scroll_cb, ...
        'KeyPressFcn', @key_cb);

    % --- Layout Configuration (3 Columns) ---
    % 1. Main Image Axis (Left ~50%)
    gui_data.pp = panel();
    gui_data.pp.pack('h', {0.7 0.3});
    gui_data.pp(2).pack('h', 2);

    gui_data.ax_pos = gui_data.pp(2,1).select();
    gui_data.ax_ang = gui_data.pp(2,2).select();
    gui_data.ax     = gui_data.pp(1).select();

    gui_data.pp.margin = [1 18 2 10];
    gui_data.pp(2).marginleft = 25;


    % gui_data.ax = axes('Parent', gui_fig, 'Position', [0.03 0.10 0.50 0.85]);
    axis(gui_data.ax, 'off', 'image');
    colormap(gui_data.ax, gray(255));
    xlim(gui_data.ax, [0.5 gui_data.Cols+0.5]);
    ylim(gui_data.ax, [0.5 gui_data.Rows+0.5]);
    
    % 2. Position Plot Axis (Middle ~20%)
    % Shared Y axis labels shown here
    % gui_data.ax_pos = axes('Parent', gui_fig, 'Position', [0.58 0.10 0.18 0.80]);
    hold(gui_data.ax_pos, 'on');
    grid(gui_data.ax_pos, 'on');
    xlabel(gui_data.ax_pos, 'Position (px)');
    ylabel(gui_data.ax_pos, 'Slice number');
    set(gui_data.ax_pos, 'YDir', 'reverse'); % Slice 1 at Top
    maxpos = max(gui_data.Rows, gui_data.Cols);
    xlim(gui_data.ax_pos, [0 maxpos]);
    xticks(gui_data.ax_pos,[0 round(maxpos/2) maxpos])
    ylim(gui_data.ax_pos, [1 gui_data.Nslices]);
 
    
    % 3. Angle Plot Axis (Right ~20%)
    % Shared Y axis, so we hide the labels here
    % gui_data.ax_ang = axes('Parent', gui_fig, 'Position', [0.80 0.10 0.18 0.80]);
    hold(gui_data.ax_ang, 'on');
    grid(gui_data.ax_ang, 'on');
    xlabel(gui_data.ax_ang, 'Angle (rad)');
    set(gui_data.ax_ang, 'YDir', 'reverse'); % Slice 1 at Top
    set(gui_data.ax_ang, 'YTickLabel', []);  % Hide Y labels
    xlim(gui_data.ax_ang, [-pi-0.1 pi+0.1]);
    ylim(gui_data.ax_ang, [1 gui_data.Nslices]);
    
    % % Link Y-axes so they stay in sync if code were to change zoom (optional but good practice)
    % linkaxes([gui_data.ax_pos, gui_data.ax_ang], 'y');

    % --- Main Image Objects ---
    gui_data.h_im = imagesc(zeros(100), 'Parent', gui_data.ax);
    set(gui_data.h_im, 'ButtonDownFcn', @mouse_cb);
    
    % User Inputs (Dots)
    gui_data.h_user_cen = line(gui_data.ax, nan, nan, 'Color','b', 'Marker','.', ...
        'MarkerSize', 25, 'PickableParts','none'); 
    gui_data.h_user_ant = line(gui_data.ax, nan, nan, 'Color','g', 'Marker', '.',...
        'MarkerSize', 25, 'PickableParts','none'); 
    gui_data.h_user_pos = line(gui_data.ax, nan, nan, 'Color','r','Marker','.', ...
        'MarkerSize', 25, 'PickableParts','none'); 
    
    % Optimization Predictions (Image Overlay)
    gui_data.h_fit_cen  = line(gui_data.ax, nan, nan, 'Color', 'c','Marker','+', ...
        'MarkerSize', 10, 'LineWidth', 2, 'PickableParts','none');
    gui_data.h_fit_line = line(gui_data.ax, nan, nan, 'Color','y','LineStyle','-', ...
        'LineWidth', 1.5, 'PickableParts','none');
    gui_data.h_fit_ant  = line(gui_data.ax, nan, nan, 'Color', 'g','Marker','s', ...
        'MarkerSize', 6, 'LineWidth', 1, 'PickableParts','none');
    gui_data.h_fit_pos  = line(gui_data.ax, nan, nan, 'Color', 'r','Marker','s', ...
        'MarkerSize', 6, 'LineWidth', 1, 'PickableParts','none');

    % --- Side Panel Objects ---
    slice_vec = 1:gui_data.Nslices;
    
    % Plot 1: Position (X=Pos, Y=Slice)
    % Fitted Lines
    gui_data.h_plot_x_line = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Fit X');
    gui_data.h_plot_y_line = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'm-', 'LineWidth', 1.2, 'DisplayName', 'Fit Y');
    % Raw User Observations (Dots)
    gui_data.h_plot_x_obs  = plot(gui_data.ax_pos, nan, nan, 'bo', 'MarkerSize', 5, 'DisplayName', 'User X', 'PickableParts','none');
    gui_data.h_plot_y_obs  = plot(gui_data.ax_pos, nan, nan, 'mo', 'MarkerSize', 5, 'DisplayName', 'User Y', 'PickableParts','none');
    % Cursor Line (Horizontal)
    gui_data.h_cursor_pos  = plot(gui_data.ax_pos, [-1e5 1e5], [1 1], 'k--', 'LineWidth', 1, 'PickableParts','none');
    legend(gui_data.ax_pos, [gui_data.h_plot_x_line gui_data.h_plot_y_line], 'Location', 'northeast');

    % Plot 2: Angle (X=Ang, Y=Slice)
    % Fitted Line
    gui_data.h_plot_th_line = plot(gui_data.ax_ang, nan(size(slice_vec)), slice_vec, 'Color', [0 0.6 0], 'LineWidth', 1.2, 'DisplayName', 'Fit \theta');
    % Raw User Observations (Dots)
    gui_data.h_plot_th_obs  = plot(gui_data.ax_ang, nan, nan, 'o', 'Color', [0 0.6 0], 'MarkerSize', 5, 'DisplayName', 'User \theta', 'PickableParts','none');
    % Cursor Line (Horizontal)
    gui_data.h_cursor_ang   = plot(gui_data.ax_ang, [-10 10], [1 1], 'k--', 'LineWidth', 1, 'PickableParts','none');
    
    % --- Link Navigation Callbacks ---
    % Attach click callback to axes and lines
    set(gui_data.ax_pos, 'ButtonDownFcn', @plot_click_cb);
    set(gui_data.h_plot_x_line, 'ButtonDownFcn', @plot_click_cb);
    set(gui_data.h_plot_y_line, 'ButtonDownFcn', @plot_click_cb);
    
    set(gui_data.ax_ang, 'ButtonDownFcn', @plot_click_cb);
    set(gui_data.h_plot_th_line, 'ButtonDownFcn', @plot_click_cb);
    
    guidata(gui_fig, gui_data);
    update_view(gui_fig);
end

% -------------------------------------------------------------------------
%   INTERACTION CALLBACKS
% -------------------------------------------------------------------------
function mouse_cb(src, ~)
    fig = ancestor(src, 'figure');
    gui_data = guidata(fig);
    
    pt = get(gui_data.ax, 'CurrentPoint');
    x = pt(1,1); y = pt(1,2);
    s = gui_data.curr_slice;
    
    click_type = get(fig, 'SelectionType'); 
    
    if strcmp(click_type, 'normal')
        % Left -> ANTERIOR (Green)
        gui_data.user.ant(s, :) = [x, y];
        
    elseif strcmp(click_type, 'alt')
        % Right -> POSTERIOR (Red)
        gui_data.user.pos(s, :) = [x, y];
        
    elseif strcmp(click_type, 'extend')
        % Middle -> CENTER CANAL (Blue)
        gui_data.user.cen(s, :) = [x, y];
    end
    
    gui_data = run_optimizer(gui_data);
    guidata(fig, gui_data);
    update_view(fig);
end

function plot_click_cb(src, ~)
    % Handle clicks on the side plots to navigate slices
    % NOTE: Plots are vertical (Y axis = slice index)
    fig = ancestor(src, 'figure');
    gui_data = guidata(fig);
    
    if isa(src, 'matlab.graphics.axis.Axes')
        ax_clicked = src;
    else
        ax_clicked = get(src, 'Parent');
    end
    
    pt = get(ax_clicked, 'CurrentPoint');
    % pt(1,2) corresponds to the Y-axis value, which is now our slice index
    new_slice_idx = round(pt(1,2));
    
    change_slice(fig, new_slice_idx);
end

function key_cb(fig, evt)
    gui_data = guidata(fig);
    s = gui_data.curr_slice;
    
    switch evt.Key
        case 'leftarrow',  change_slice(fig, s - 1);
        case 'rightarrow', change_slice(fig, s + 1);
        
        case 'l'
            % Change Regularization Parameters
            prompt = {'Positional Lambda:', 'Angular Lambda:'};
            dlgtitle = 'Regularization Settings';
            dims = [1 35];
            definput = {num2str(gui_data.lambda_pos), num2str(gui_data.lambda_ang)};
            
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                new_pos = str2double(answer{1});
                new_ang = str2double(answer{2});
                
                if ~isnan(new_pos) && ~isnan(new_ang)
                    gui_data.lambda_pos = new_pos;
                    gui_data.lambda_ang = new_ang;
                    % Re-run solver with new values
                    gui_data = run_optimizer(gui_data);
                    guidata(fig, gui_data);
                    update_view(fig);
                end
            end

        case 'p'
            % Toggle Prediction Visibility
            gui_data.show_pred = ~gui_data.show_pred;
            guidata(fig, gui_data);
            update_view(fig);
            
        case 'c'
            gui_data.user.cen(s,:) = [nan nan];
            gui_data.user.ant(s,:) = [nan nan];
            gui_data.user.pos(s,:) = [nan nan];
            gui_data = run_optimizer(gui_data);
            guidata(fig, gui_data);
            update_view(fig);
            
        case 'x' 
            pt = get(gui_data.ax, 'CurrentPoint'); 
            pt = pt(1,1:2);
            d_cen = norm(gui_data.user.cen(s,:) - pt);
            d_ant = norm(gui_data.user.ant(s,:) - pt);
            d_pos = norm(gui_data.user.pos(s,:) - pt);
            
            if isnan(d_cen), d_cen = inf; end
            if isnan(d_ant), d_ant = inf; end
            if isnan(d_pos), d_pos = inf; end
            
            [mval, midx] = min([d_cen, d_ant, d_pos]);
            if mval < 50 
                if midx==1, gui_data.user.cen(s,:) = [nan nan]; end
                if midx==2, gui_data.user.ant(s,:) = [nan nan]; end
                if midx==3, gui_data.user.pos(s,:) = [nan nan]; end
                gui_data = run_optimizer(gui_data);
                guidata(fig, gui_data);
                update_view(fig);
            end
        case 'space'
            has_data = ~isnan(gui_data.user.cen(:,1)) | ...
                       ~isnan(gui_data.user.ant(:,1)) | ...
                       ~isnan(gui_data.user.pos(:,1));
            
            if ~any(has_data)
                change_slice(fig, round(gui_data.Nslices/2));
                return;
            end
            
            dist_map = bwdist(has_data);
            dist_map(s) = 0; 
            [~, idx] = max(dist_map);
            change_slice(fig, idx);
            
        case 's'
            save_data(gui_data);
    end
end

function scroll_cb(fig, evt)
    gui_data = guidata(fig);
    change_slice(fig, gui_data.curr_slice + evt.VerticalScrollCount);
end

function change_slice(fig, new_idx)
    gui_data = guidata(fig);
    gui_data.curr_slice = max(1, min(gui_data.Nslices, new_idx));
    guidata(fig, gui_data);
    update_view(fig);
end

% -------------------------------------------------------------------------
%   VISUALIZATION
% -------------------------------------------------------------------------
function update_view(fig)
    gui_data = guidata(fig);
    s = gui_data.curr_slice;
    slice_vec = 1:gui_data.Nslices;
    
    % 1. Image
    set(gui_data.h_im, 'CData', gui_data.disp_vol(:,:,s));
    
    % 2. User Dots on Image
    set(gui_data.h_user_cen, 'XData', gui_data.user.cen(s,1), 'YData', gui_data.user.cen(s,2));
    set(gui_data.h_user_ant, 'XData', gui_data.user.ant(s,1), 'YData', gui_data.user.ant(s,2));
    set(gui_data.h_user_pos, 'XData', gui_data.user.pos(s,1), 'YData', gui_data.user.pos(s,2));
    
    % 3. Model Predictions on Image
    if gui_data.show_pred && ~isnan(gui_data.fit.x(s))
        cx = gui_data.fit.x(s);
        cy = gui_data.fit.y(s);
        th = gui_data.fit.theta(s);
        r  = gui_data.fit.rad(s);
        if isnan(r), r = 10; end 
        
        dx = r * cos(th);
        dy = r * sin(th);
        
        set(gui_data.h_fit_cen, 'XData', cx, 'YData', cy, 'Visible', 'on');
        set(gui_data.h_fit_line, 'XData', [cx-dx, cx+dx], 'YData', [cy-dy, cy+dy], 'Visible', 'on');
        set(gui_data.h_fit_ant, 'XData', cx+dx, 'YData', cy+dy, 'Visible', 'on');
        set(gui_data.h_fit_pos, 'XData', cx-dx, 'YData', cy-dy, 'Visible', 'on');
    else
        set(gui_data.h_fit_cen, 'Visible', 'off');
        set(gui_data.h_fit_line, 'Visible', 'off');
        set(gui_data.h_fit_ant, 'Visible', 'off');
        set(gui_data.h_fit_pos, 'Visible', 'off');
    end
    
    % 4. Update Side Panels (Plots)
    % Note: XData = Value, YData = Slice Index (Vertical orientation)
    
    % Position Plot
    set(gui_data.h_plot_x_line, 'XData', gui_data.fit.x, 'YData', slice_vec);
    set(gui_data.h_plot_y_line, 'XData', gui_data.fit.y, 'YData', slice_vec);
    set(gui_data.h_plot_x_obs,  'XData', gui_data.obs.x, 'YData', slice_vec);
    set(gui_data.h_plot_y_obs,  'XData', gui_data.obs.y, 'YData', slice_vec);
    set(gui_data.h_cursor_pos,  'YData', [s s]); % Horizontal line at slice s
    
    % Angle Plot
    set(gui_data.h_plot_th_line, 'XData', gui_data.fit.theta, 'YData', slice_vec);
    set(gui_data.h_plot_th_obs,  'XData', gui_data.obs.th,    'YData', slice_vec);
    set(gui_data.h_cursor_ang,   'YData', [s s]);
        
    % 5. Title
    n_cen = sum(~isnan(gui_data.user.cen(:,1)));
    n_ant = sum(~isnan(gui_data.user.ant(:,1)));
    n_pos = sum(~isnan(gui_data.user.pos(:,1)));
    
    status_str = '';
    if ~gui_data.show_pred, status_str = ' [PREDICTION HIDDEN]'; end
    
    key_str = "[LClick] Ant | [RClick] Pos | [MidClick] Cen | [c] Clear | [l] Set Reg | [p] Show Fit | [s] Save";
    tstr = sprintf('Slice: %d | Clicks: C=%d A=%d P=%d | Regul.: Pos=%.0f Ang=%.0f%s\n%s', ...
        s, n_cen, n_ant, n_pos, gui_data.lambda_pos, gui_data.lambda_ang, status_str, key_str);
    title(gui_data.ax, tstr, 'Color', 'k', 'FontSize', 11, 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
%   CORE OPTIMIZATION SOLVER
% -------------------------------------------------------------------------
function gui_data = run_optimizer(gui_data)
    N = gui_data.Nslices;
    
    % Prepare observations
    obs_x = nan(N, 1);
    obs_y = nan(N, 1);
    obs_th = nan(N, 1);
    obs_rad = nan(N, 1);
    
    for z = 1:N
        c = gui_data.user.cen(z,:);
        a = gui_data.user.ant(z,:);
        p = gui_data.user.pos(z,:);
        
        has_c = ~isnan(c(1));
        has_a = ~isnan(a(1));
        has_p = ~isnan(p(1));
        
        % A. Center
        if has_c
            obs_x(z) = c(1); obs_y(z) = c(2);
        elseif has_a && has_p
            obs_x(z) = (a(1) + p(1))/2; obs_y(z) = (a(2) + p(2))/2;
        end
        
        % B. Angle & Radius
        vec = [nan nan]; curr_rad = nan;
        if has_a && has_p
            vec = a - p; curr_rad = norm(vec)/2;
        elseif has_c && has_a
            vec = a - c; curr_rad = norm(vec);
        elseif has_c && has_p
            vec = c - p; curr_rad = norm(vec);
        end
        
        if ~isnan(vec(1))
            obs_th(z) = atan2(vec(2), vec(1));
            obs_rad(z) = curr_rad;
        end
    end
    
    % Store raw observations in gui_data for plotting
    gui_data.obs.x = obs_x;
    gui_data.obs.y = obs_y;
    gui_data.obs.th = obs_th;
    
    % Solver (X, Y)
    idx_x = find(~isnan(obs_x));
    if length(idx_x) > 1
        gui_data.fit.x = full(solve_spline(N, idx_x, obs_x(idx_x), gui_data.lambda_pos));
        gui_data.fit.y = full(solve_spline(N, idx_x, obs_y(idx_x), gui_data.lambda_pos));
    else
        gui_data.fit.x = nan(N,1); gui_data.fit.y = nan(N,1);
    end
    
    % Solver (Theta)
    idx_th = find(~isnan(obs_th));
    if length(idx_th) > 1
        raw = obs_th(idx_th);
        unwrapped = unwrap_sparse(raw);
        fit_unwrapped = solve_spline(N, idx_th, unwrapped, gui_data.lambda_ang);
        gui_data.fit.theta = full(angle(exp(1i * fit_unwrapped))); 
    else
        gui_data.fit.theta = nan(N,1);
    end
    
    % Solver (Radius)
    idx_r = find(~isnan(obs_rad));
    if length(idx_r) > 1
        gui_data.fit.rad = solve_spline(N, idx_r, obs_rad(idx_r), gui_data.lambda_pos);
    else
        gui_data.fit.rad = repmat(mean(obs_rad, 'omitnan'), N, 1);
    end
end

function x_full = solve_spline(N, obs_idx, obs_val, lambda)
    e = ones(N, 1);
    D2 = spdiags([e -2*e e], -1:1, N, N);
    D2(1,:) = 0; D2(1,1:2) = [1 -1];
    D2(N,:) = 0; D2(N,N-1:N) = [-1 1];
    L = D2' * D2;
    W = sparse(obs_idx, obs_idx, ones(length(obs_idx),1), N, N);
    b = sparse(N, 1); b(obs_idx) = obs_val;
    x_full = (W + lambda * L) \ b;
end

function th_unwrapped = unwrap_sparse(th)
    th_unwrapped = th;
    for k = 2:length(th)
        diff = th_unwrapped(k) - th_unwrapped(k-1);
        shift = round(diff / (2*pi));
        th_unwrapped(k) = th_unwrapped(k) - shift * 2 * pi;
    end
end

function save_data(gui_data)
    align_out.user_cen = gui_data.user.cen;
    align_out.user_ant = gui_data.user.ant;
    align_out.user_pos = gui_data.user.pos;
    align_out.fit_x = gui_data.fit.x;
    align_out.fit_y = gui_data.fit.y;
    align_out.fit_theta = gui_data.fit.theta;
    
    align_out.lambda_pos = gui_data.lambda_pos;
    align_out.lambda_ang = gui_data.lambda_ang;
    
    fname = fullfile(gui_data.save_path, 'spinal_alignment_opt.mat');
    save(fname, 'align_out');
    
    t = get(gui_data.ax, 'Title');
    old_t = t.String;
    if ischar(old_t)
        new_t = [old_t, ' [SAVED]'];
    else
        new_t = old_t;
        new_t{end} = [new_t{end}, ' [SAVED]'];
    end
    
    title(gui_data.ax, new_t, 'Color', 'g');
    pause(0.2);
    update_view(ancestor(gui_data.ax, 'figure'));
end