function spinal_cord_aligner(opts)
%SPINAL_CORD_ALIGNER Refine the centre line and twist of a spinal cord sample.
%
%   SPINAL_CORD_ALIGNER(OPTS) opens the GUI that decides how every section of
%   the cord is translated and rotated before the atlas is fitted to it. OPTS is
%   the options struct written by PREPARECORDSAMPLEFORREGISTRATION, which is
%   also the second half of the initialisation: it already carries a *prediction*
%   of the centre and the dorsoventral angle of every section, computed from the
%   Otsu segmentation and the direction of intensity variance within each plane
%   (see AUTOCORDCENTERLINE).
%
%   So this GUI is a correction tool, not an annotation tool. It starts with a
%   complete centre line and fits one smoothing spline to a mixed set of
%   observations: the points the user clicks, and the automatic prediction on
%   every slice the user has not spoken for.
%
%   A click speaks for a stretch of cord, not just for its own slice. On the
%   clicked slice and OPTS.userexclusionradius slices either side (15 by
%   default) the prediction is dropped from the fit completely, so the curve
%   there follows the user's points and nothing else. That is what makes a
%   correction stick: the prediction is an observation on nearly every slice, so
%   a click competing against it directly would hardly move the curve. Away from
%   any click the prediction is followed as before, and clicking nothing at all
%   is a valid answer when the prediction is already right.
%
%   The dashed prediction in the side plots is drawn only where it is still
%   being used, so the stretches the user has taken over show up as gaps.
%
%   The one thing the prediction cannot know is which end of the intensity axis
%   is anterior - that depends on the stain, not on the image. If the whole cord
%   comes out rotated by 180 degrees, press 'f' once.
%
%   CONTROLS
%       Left click   : set ANTERIOR/FRONT (green)
%       Right click  : set POSTERIOR/BACK (red)
%       Middle click : set CENTER CANAL (blue)
%       'f'          : flip the predicted angle by 180 degrees, whole cord
%       'l'          : set the regularisation and how far a click reaches
%       'p'          : toggle prediction visibility
%       'c'          : clear all points on the current slice
%       'x'          : delete the point nearest the cursor
%       Space        : jump to the largest gap between corrections
%       's'          : save
%       Click on the side plots to jump to a slice
%
%   The result is saved as 'spinal_alignment_opt.mat' in OPTS.savepath and read
%   back by INITIALIZECORDREGISTRATION.
%
%   See also PREPARECORDSAMPLEFORREGISTRATION, AUTOCORDCENTERLINE,
%   INITIALIZECORDREGISTRATION, COMPUTESTRAIGHTENINGTRANSFORMS.

    % -- 1. Setup & Input Check --
    if nargin < 1 || ~isstruct(opts)
        error('spinal_cord_aligner:badInput', ...
            'Input must be the options struct from prepareCordSampleForRegistration.');
    end

    gui_data = struct();
    gui_data.save_path = getOr(opts, {'savepath', 'lsfolder'}, '');
    if isempty(gui_data.save_path)
        error('spinal_cord_aligner:noSavePath', ...
            'opts must carry a savepath (or a legacy lsfolder) to save into.');
    end
    if ~isfolder(gui_data.save_path)
        mkdir(gui_data.save_path);
    end

    % Volume: written to disk by prepareCordSampleForRegistration, but an
    % in-memory volume is still accepted for one-off use.
    if isfield(opts, 'regvol') && ~isempty(opts.regvol)
        regvol = opts.regvol;
    elseif isfield(opts, 'cordvolpath') && exist(opts.cordvolpath, 'file') == 2
        fprintf('Reading cord volume %s\n', opts.cordvolpath);
        regvol = readDownStack(opts.cordvolpath);
    else
        error('spinal_cord_aligner:noVolume', ...
            ['No cord volume found. Run prepareCordSampleForRegistration first, ' ...
             'so that opts.cordvolpath points at the cropped registration volume.']);
    end

    [gui_data.Rows, gui_data.Cols, gui_data.Nslices] = size(regvol);

    % Contrast Normalization
    rng(1);
    samp_idx = randperm(numel(regvol), min(numel(regvol), 2e4));
    v_samp   = single(regvol(samp_idx));
    vmax     = quantile(v_samp, 0.999);
    vmin     = quantile(v_samp, 0.001);
    v_float  = (single(regvol) - vmin) / max(vmax - vmin, eps);
    v_float  = max(0, min(1, v_float));
    gui_data.disp_vol = uint8(255 * v_float);
    clear regvol v_float;

   % -- Automatic predictions (the starting point of every fit) --
    gui_data.auto = parseAutoPrediction(opts, gui_data.Nslices);
    gui_data.theta_flip = false;

   % -- Data Structures --
    gui_data.user.cen = nan(gui_data.Nslices, 2);
    gui_data.user.ant = nan(gui_data.Nslices, 2);
    gui_data.user.pos = nan(gui_data.Nslices, 2);

    % Derived Observations (for plotting raw user points)
    gui_data.obs.x  = nan(gui_data.Nslices, 1);
    gui_data.obs.y  = nan(gui_data.Nslices, 1);
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

    % How far a user point reaches, in slices: the prediction is dropped from
    % the fit on the clicked slice and this many slices either side of it. It is
    % wide because the regularisation is stiff - a click needs a stretch of cord
    % to itself before it can bend the curve at all.
    gui_data.excl_radius  = getOr(opts, 'userexclusionradius', 100);
    gui_data.autoused.pos = true(gui_data.Nslices, 1);
    gui_data.autoused.ang = true(gui_data.Nslices, 1);

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
                    if isfield(d, 'theta_flip')
                        gui_data.theta_flip = d.theta_flip;
                    end
                    if isfield(d, 'excl_radius')
                        gui_data.excl_radius = d.excl_radius;
                    end

                    fprintf('Data loaded (Reg: Pos=%.0f, Ang=%.0f, Reach=%d slices).\n', ...
                        gui_data.lambda_pos, gui_data.lambda_ang, gui_data.excl_radius);
                end
            catch
                warning('Load failed.');
            end
        end
    end
    gui_data = run_optimizer(gui_data);

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

    axis(gui_data.ax, 'off', 'image');
    colormap(gui_data.ax, gray(255));
    xlim(gui_data.ax, [0.5 gui_data.Cols+0.5]);
    ylim(gui_data.ax, [0.5 gui_data.Rows+0.5]);

    % 2. Position Plot Axis (Middle ~20%)
    % Shared Y axis labels shown here
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
    hold(gui_data.ax_ang, 'on');
    grid(gui_data.ax_ang, 'on');
    xlabel(gui_data.ax_ang, 'Angle (rad)');
    set(gui_data.ax_ang, 'YDir', 'reverse'); % Slice 1 at Top
    set(gui_data.ax_ang, 'YTickLabel', []);  % Hide Y labels
    xlim(gui_data.ax_ang, [-pi-0.1 pi+0.1]);
    ylim(gui_data.ax_ang, [1 gui_data.Nslices]);

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
    % Automatic predictions (dashed), fitted lines (solid), user clicks (dots)
    gui_data.h_plot_x_auto = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'b--', 'LineWidth', 0.8, 'DisplayName', 'Auto X');
    gui_data.h_plot_y_auto = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'm--', 'LineWidth', 0.8, 'DisplayName', 'Auto Y');
    gui_data.h_plot_x_line = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Fit X');
    gui_data.h_plot_y_line = plot(gui_data.ax_pos, nan(size(slice_vec)), slice_vec, 'm-', 'LineWidth', 1.2, 'DisplayName', 'Fit Y');
    gui_data.h_plot_x_obs  = plot(gui_data.ax_pos, nan, nan, 'bo', 'MarkerSize', 5, 'DisplayName', 'User X', 'PickableParts','none');
    gui_data.h_plot_y_obs  = plot(gui_data.ax_pos, nan, nan, 'mo', 'MarkerSize', 5, 'DisplayName', 'User Y', 'PickableParts','none');
    % Cursor Line (Horizontal)
    gui_data.h_cursor_pos  = plot(gui_data.ax_pos, [-1e5 1e5], [1 1], 'k--', 'LineWidth', 1, 'PickableParts','none');
    legend(gui_data.ax_pos, [gui_data.h_plot_x_line gui_data.h_plot_y_line], 'Location', 'northeast');

    % Plot 2: Angle (X=Ang, Y=Slice)
    gui_data.h_plot_th_auto = plot(gui_data.ax_ang, nan(size(slice_vec)), slice_vec, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'Auto \theta');
    gui_data.h_plot_th_line = plot(gui_data.ax_ang, nan(size(slice_vec)), slice_vec, 'Color', [0 0.6 0], 'LineWidth', 1.2, 'DisplayName', 'Fit \theta');
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
            prompt = {'Positional Lambda:', 'Angular Lambda:', ...
                'Reach of a user point (slices, prediction dropped within):'};
            dlgtitle = 'Regularization Settings';
            dims = [1 55];
            definput = {num2str(gui_data.lambda_pos), num2str(gui_data.lambda_ang), ...
                num2str(gui_data.excl_radius)};

            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                new_pos = str2double(answer{1});
                new_ang = str2double(answer{2});
                new_rad = str2double(answer{3});

                if ~isnan(new_pos) && ~isnan(new_ang) && ~isnan(new_rad) && new_rad >= 0
                    gui_data.lambda_pos  = new_pos;
                    gui_data.lambda_ang  = new_ang;
                    gui_data.excl_radius = round(new_rad);
                    % Re-run solver with new values
                    gui_data = run_optimizer(gui_data);
                    guidata(fig, gui_data);
                    update_view(fig);
                end
            end

        case 'f'
            % Flip the automatic angle prediction by 180 degrees for the whole
            % cord. The image alone cannot say which end of the intensity axis
            % is anterior, so this is the one global choice left to the user.
            if isempty(gui_data.auto.theta)
                return;
            end
            gui_data.theta_flip = ~gui_data.theta_flip;
            gui_data = run_optimizer(gui_data);
            guidata(fig, gui_data);
            update_view(fig);

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
    auto_th   = autoTheta(gui_data);

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
        if isnan(r) || r <= 0, r = 10; end

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
    % the automatic series is drawn only where it is actually used as an
    % observation, so the stretches the user has taken over are visible as gaps
    set(gui_data.h_plot_x_auto, 'XData', usedSeries(gui_data.auto.x, gui_data.autoused.pos, gui_data.Nslices), 'YData', slice_vec);
    set(gui_data.h_plot_y_auto, 'XData', usedSeries(gui_data.auto.y, gui_data.autoused.pos, gui_data.Nslices), 'YData', slice_vec);
    set(gui_data.h_plot_x_line, 'XData', gui_data.fit.x, 'YData', slice_vec);
    set(gui_data.h_plot_y_line, 'XData', gui_data.fit.y, 'YData', slice_vec);
    set(gui_data.h_plot_x_obs,  'XData', gui_data.obs.x, 'YData', slice_vec);
    set(gui_data.h_plot_y_obs,  'XData', gui_data.obs.y, 'YData', slice_vec);
    set(gui_data.h_cursor_pos,  'YData', [s s]); % Horizontal line at slice s

    % Angle Plot
    set(gui_data.h_plot_th_auto, 'XData', usedSeries(auto_th, gui_data.autoused.ang, gui_data.Nslices), 'YData', slice_vec);
    set(gui_data.h_plot_th_line, 'XData', gui_data.fit.theta, 'YData', slice_vec);
    set(gui_data.h_plot_th_obs,  'XData', gui_data.obs.th,    'YData', slice_vec);
    set(gui_data.h_cursor_ang,   'YData', [s s]);

    % 5. Title
    n_cen = sum(~isnan(gui_data.user.cen(:,1)));
    n_ant = sum(~isnan(gui_data.user.ant(:,1)));
    n_pos = sum(~isnan(gui_data.user.pos(:,1)));

    status_str = '';
    if ~gui_data.show_pred, status_str = [status_str ' [PREDICTION HIDDEN]']; end
    if isempty(gui_data.auto.x)
        status_str = [status_str ' [NO AUTO PREDICTION]'];
    elseif gui_data.theta_flip
        status_str = [status_str ' [ANGLE FLIPPED]'];
    end

    % slices where the prediction has been dropped in favour of the user
    n_owned = nnz(~gui_data.autoused.pos);

    key_str = "[LClick] Ant | [RClick] Pos | [MidClick] Cen | [f] Flip 180 | [c] Clear | [l] Set Reg | [p] Show Fit | [s] Save";
    tstr = sprintf(['Slice: %d | Clicks: C=%d A=%d P=%d | Regul.: Pos=%.0f Ang=%.0f | ' ...
        'Reach: %d slices (prediction dropped on %d/%d)%s\n%s'], ...
        s, n_cen, n_ant, n_pos, gui_data.lambda_pos, gui_data.lambda_ang, ...
        gui_data.excl_radius, n_owned, gui_data.Nslices, status_str, key_str);
    title(gui_data.ax, tstr, 'Color', 'k', 'FontSize', 11, 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
%   CORE OPTIMIZATION SOLVER
% -------------------------------------------------------------------------
function gui_data = run_optimizer(gui_data)
%RUN_OPTIMIZER Fit one spline through the user's points and the predictions.
%   There is a single smoothing spline per quantity, and it is fitted to a mixed
%   set of observations: the points the user clicked, plus the automatic
%   prediction of AUTOCORDCENTERLINE at every slice the user has *not* spoken
%   for.
%
%   A click speaks for its own slice and for gui_data.excl_radius slices on
%   either side: inside that window the automatic prediction is dropped from the
%   fit entirely, so the curve there is decided by the user's points and the
%   smoothness penalty alone. Outside it the prediction is an observation like
%   any other, which is why a cord nobody has clicked on comes out as the
%   prediction itself.
%
%   That window is what makes a correction stick. The prediction is an
%   observation on nearly every slice, so a click competing against it directly
%   would barely move the curve; removing its neighbours first is what lets a
%   handful of clicks take over a stretch of cord. Widen the window with 'l' if
%   a correction is still being pulled back by the prediction next to it.
%
%   The exclusion is per quantity, not per slice: a middle click gives a centre
%   but says nothing about the angle, so it suppresses the predicted centres
%   around it and leaves the predicted angles alone.
%
%   Angles are fitted as vectors - the cosine and the sine are smoothed and the
%   angle read back off the result - so nothing anywhere has to choose a branch
%   for an observation. Two clicks a degree apart are a degree apart no matter
%   where the +/-pi seam falls or how far the prediction is from them.
%
%   One thing to expect from a second-difference spline: between two clicks that
%   both correct in the same direction, the curve overshoots them a little,
%   because it still has to swing back to the prediction just outside the
%   window. It shrinks as the reach or the regularisation is raised (both under
%   'l'), and it does not appear around a single click - at the defaults a lone
%   click lands about 92% of the way to where it was put, decaying to nothing by
%   the edge of its window.

    N       = gui_data.Nslices;
    auto_th = autoTheta(gui_data);
    hasauto = ~isempty(gui_data.auto.x);
    radius  = gui_data.excl_radius;

    % Prepare observations
    obs_x   = nan(N, 1);
    obs_y   = nan(N, 1);
    obs_th  = nan(N, 1);
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
            obs_th(z)  = atan2(vec(2), vec(1));
            obs_rad(z) = curr_rad;
        end
    end

    % Store raw observations in gui_data for plotting
    gui_data.obs.x  = obs_x;
    gui_data.obs.y  = obs_y;
    gui_data.obs.th = obs_th;

    % Solver (X, Y): user centres, plus predicted centres out of their reach
    idx_u = find(~isnan(obs_x));
    [idx_x, val_x, keep_x] = mixCordObservations(N, idx_u, obs_x(idx_u), gui_data.auto.x, radius);
    [~,     val_y]         = mixCordObservations(N, idx_u, obs_y(idx_u), gui_data.auto.y, radius);
    gui_data.fit.x = fitCordSeries(N, idx_x, val_x, gui_data.lambda_pos);
    gui_data.fit.y = fitCordSeries(N, idx_x, val_y, gui_data.lambda_pos);
    gui_data.autoused.pos = keep_x;

    % Solver (Theta): fitted as a vector, never as an angle.
    %
    % Angles have no ordering to smooth along, so anything that fits them
    % directly has to first pick a branch for each observation, and picking it
    % relative to the prediction breaks down exactly where the two disagree by
    % about half a turn: a hair's difference between two clicks then puts them a
    % full turn apart and the fit swings wildly between them. Smoothing the
    % cosine and the sine instead has no branch to pick and no seam to cross -
    % two clicks a degree apart are always a degree apart - and the angle comes
    % back out of the fitted vector at the end.
    idx_u = find(~isnan(obs_th));
    if hasauto
        autoc = cos(auto_th);
        autos = sin(auto_th);
    else
        autoc = [];
        autos = [];
    end
    [idx_th, val_c, keep_th] = mixCordObservations(N, idx_u, cos(obs_th(idx_u)), autoc, radius);
    [~,      val_s]          = mixCordObservations(N, idx_u, sin(obs_th(idx_u)), autos, radius);

    fitc = fitCordSeries(N, idx_th, val_c, gui_data.lambda_ang);
    fits = fitCordSeries(N, idx_th, val_s, gui_data.lambda_ang);
    gui_data.fit.theta    = atan2(fits, fitc);
    gui_data.autoused.ang = keep_th;

    % Solver (Radius): display only
    idx_u = find(~isnan(obs_rad));
    [idx_r, val_r] = mixCordObservations(N, idx_u, obs_rad(idx_u), gui_data.auto.rad, radius);
    gui_data.fit.rad = fitCordSeries(N, idx_r, val_r, gui_data.lambda_pos);
    if all(isnan(gui_data.fit.rad))
        gui_data.fit.rad = repmat(mean(obs_rad, 'omitnan'), N, 1);
    end
end


function th = autoTheta(gui_data)
%AUTOTHETA Predicted angle for display, wrapped, with the global flip applied.
    th = gui_data.auto.theta;
    if ~isempty(th) && gui_data.theta_flip
        th = wrapToPiLocal(th + pi);
    end
end

function v = usedSeries(v, used, N)
%USEDSERIES The predicted series, blanked wherever it is not used in the fit.
    if isempty(v)
        v = nan(N, 1);
        return
    end
    v(~used) = NaN;
end

function auto = parseAutoPrediction(opts, Nslices)
%PARSEAUTOPREDICTION Pull the automatic centre line out of the options struct.
    auto = struct('x', [], 'y', [], 'theta', [], 'thetaunwr', [], 'rad', []);
    cordauto = getOr(opts, 'cordauto', []);
    if isempty(cordauto)
        warning('spinal_cord_aligner:noPrediction', ...
            ['No automatic centre line in the options struct, so every slice has ' ...
             'to be defined by hand. Run prepareCordSampleForRegistration to get ' ...
             'predictions to correct instead.']);
        return
    end
    if numel(cordauto.cen_x) ~= Nslices
        warning('spinal_cord_aligner:predictionSizeMismatch', ...
            ['The automatic centre line covers %d slices but the volume has %d; ' ...
             'ignoring it. Re-run prepareCordSampleForRegistration.'], ...
            numel(cordauto.cen_x), Nslices);
        return
    end
    auto.x     = cordauto.cen_x(:);
    auto.y     = cordauto.cen_y(:);
    auto.theta = cordauto.theta(:);
    auto.rad   = cordauto.rad(:);
    % predictions written before the angle was tracked continuously carry no
    % continuous branch, so fall back to unwrapping the wrapped one
    auto.thetaunwr = getOr(cordauto, 'thetaunwr', unwrap(cordauto.theta(:)));
    auto.thetaunwr = auto.thetaunwr(:);
    fprintf('Starting from the automatic centre line (%d slices, net twist %2.0f deg).\n', ...
        Nslices, rad2deg(auto.thetaunwr(end) - auto.thetaunwr(1)));
end

function save_data(gui_data)
    align_out.user_cen   = gui_data.user.cen;
    align_out.user_ant   = gui_data.user.ant;
    align_out.user_pos   = gui_data.user.pos;
    align_out.fit_x      = gui_data.fit.x;
    align_out.fit_y      = gui_data.fit.y;
    align_out.fit_theta  = gui_data.fit.theta;
    align_out.fit_rad    = gui_data.fit.rad;

    % what the fit started from, so the correction can be told apart from the
    % prediction later on
    align_out.auto_x     = gui_data.auto.x;
    align_out.auto_y     = gui_data.auto.y;
    align_out.auto_theta = autoTheta(gui_data);
    align_out.theta_flip = gui_data.theta_flip;

    % which slices the fit took from the prediction and which from the user
    align_out.auto_used_pos = gui_data.autoused.pos;
    align_out.auto_used_ang = gui_data.autoused.ang;

    align_out.lambda_pos  = gui_data.lambda_pos;
    align_out.lambda_ang  = gui_data.lambda_ang;
    align_out.excl_radius = gui_data.excl_radius;

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
