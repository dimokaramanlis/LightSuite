function matchControlPoints_unified(opts)
% MATCHCONTROLPOINTS_UNIFIED Unified GUI for alignment of histology and CCF.
%   Merges functionality of matchControlPoints (Brain) and matchControlPointsSpine.
%
%   Usage:
%       matchControlPoints_unified(regopts) % Spine Mode
%       matchControlPoints_unified(opts)    % Brain Mode
%
%   Part of AP_histology toolbox

% Initialize guidata
gui_data = struct;

%==========================================================================
% 1. Mode Detection & Data Loading
%==========================================================================

if isfield(opts, 'straightvol')
    mode = 'spine';
    gui_data.mode = 'spine';
    gui_data.save_path = opts.lsfolder;
    gui_data.save_filename = 'corresponding_points.mat';
    disp('Starting Control Point Matching in SPINE Mode...');
    
elseif isfield(opts, 'downfac_reg')
    mode = 'brain';
    gui_data.mode = 'brain';
    gui_data.save_path = opts.savepath;
    gui_data.save_filename = 'atlas2histology_tform.mat';
    disp('Starting Control Point Matching in BRAIN Mode...');
else
    error('Input "opts" structure does not match known Brain or Spine formats.');
end

%==========================================================================
% 2. Data Preparation (Mode Specific)
%==========================================================================

switch mode
    case 'spine'
        % --- Spine Pipeline ---
        
        % Filter and normalize template
        gui_data.tv = medfilt3(opts.tv, [3 3 3]);
        factv = 255/single(max(gui_data.tv,[],"all"));
        gui_data.tv = uint8(single(gui_data.tv)*factv);
        gui_data.av = opts.av;
        
        Rvolume = imref3d(size(opts.straightvol));
        gui_data.Rvolume = Rvolume;
        Ratlas  = imref3d(size(gui_data.tv));
        
        % Warp initial volumes to alignment space (Atlas -> Sample)
        gui_data.tv = imwarp(gui_data.tv, Ratlas, opts.affine_atlas_to_samp, 'OutputView',Rvolume);
        gui_data.av = imwarp(gui_data.av, Ratlas, opts.affine_atlas_to_samp, 'OutputView',Rvolume);
        gui_data.Rmoving = imref3d(size(gui_data.av));
        
        % Prepare Sample Volume
        gui_data.volume = opts.straightvol;
        factv = 255/single(quantile(gui_data.volume,0.999,"all"));
        gui_data.volume = uint8(single(gui_data.volume)*factv);
        
        % Generate Slice List (Spine specific randomized selection)
        slicescheck = round(linspace(1, size(opts.straightvol, 3), 100))';
        chooselist  = [slicescheck 3*ones(size(slicescheck)) ones(size(slicescheck)) ones(size(slicescheck))];
        rng(1);
        gui_data.chooselist = chooselist(randperm(numel(slicescheck)), :);
        
        disp('Spine Data loaded and filtered.')
        
    case 'brain'
        % --- Brain Pipeline ---
        
        % Load Atlas (Standard Allen CCF Logic)
        allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
        if isempty(allen_atlas_path)
            error('No CCF atlas found (add CCF atlas to path)')
        end
        disp('Loading Allen CCF atlas...');
        
        tv_raw = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
        gui_data.tv = imresize3(tv_raw, opts.downfac_reg);

        % factv  = 255/single(max(tv_raw,[],"all"));
        factv = 255/single(max(gui_data.tv,[],"all"));
        gui_data.tv = uint8(single( gui_data.tv)*factv);
        
        av_raw = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
        gui_data.av = imresize3(av_raw, opts.downfac_reg, "Method","nearest");
        
        gui_data.Rmoving = imref3d(size(gui_data.av));
        
        % Load and Prepare Sample Volume
        volume_dir = fullfile(opts.regvolpath);
        if ~isfield(opts, 'regvolpath')
             error('No regvolpath in the options...');
        end
        volload = readDownStack(volume_dir);
        volload = permute(volload, opts.permute_sample_to_atlas);
        
        % Warp Sample (Note: Brain mode warps Sample using original_trans initially)
        % Note: We store the warped volume in gui_data.volume to match Spine logic for display
        Rvolume_orig = imref3d(size(volload));
        gui_data.volume = imwarp(volload, Rvolume_orig, opts.original_trans, 'OutputView',gui_data.Rmoving);
        
        irand = randperm(numel(gui_data.volume), min(numel(gui_data.volume), 1e4));
        factv = 255/single(quantile(gui_data.volume(irand),0.999,"all"));
        gui_data.volume = uint8(single(gui_data.volume)*factv);
        
        gui_data.Rvolume = imref3d(size(gui_data.volume));
        rng(1);
        % Generate Slice List (Brain specific function)
        gui_data.chooselist = generate_cp_list_alt(gui_data.volume);
        
        disp('Brain Data loaded and processed.');
end

%==========================================================================
% 3. Load / Initialize Control Points
%==========================================================================

load_fn = fullfile(gui_data.save_path, gui_data.save_filename);

if exist(load_fn, 'file')
    oldtform = load(load_fn);
    
    % Check for transform field name variations
    if isfield(oldtform, 'atlas2histology_tform')
        gui_data.histology_ccf_auto_alignment = oldtform.atlas2histology_tform;
    end
    
    % Handle Control Points and ensure 4-column format (x,y,z,timestamp)
    % for both modes to support the unified GUI features.
    if isfield(oldtform, 'histology_control_points')
        gui_data.histology_control_points = check_point_dims(oldtform.histology_control_points);
        gui_data.atlas_control_points     = check_point_dims(oldtform.atlas_control_points);
    else
        % Fallback if file exists but points are missing
        gui_data.histology_control_points = repmat({zeros(0,4)},size(gui_data.chooselist, 1),1);
        gui_data.atlas_control_points     = repmat({zeros(0,4)},size(gui_data.chooselist, 1),1);
    end
else
    % Initialize alignment control points (N x 4: x, y, z, timestamp)
    gui_data.histology_control_points = repmat({zeros(0,4)},size(gui_data.chooselist, 1),1);
    gui_data.atlas_control_points     = repmat({zeros(0,4)},size(gui_data.chooselist, 1),1);
end

%==========================================================================
% 4. GUI Setup (Unified: Based on Spine Version)
%==========================================================================

screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7; 
gui_width_fraction = 0.6;
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... 
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... 
    gui_width_px,gui_width_px/gui_aspect_ratio];

gui_fig = figure('KeyPressFcn',@keypress, ...
    'WindowScrollWheelFcn',@scroll_atlas_slice,...
    'Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position, ...
    'CloseRequestFcn',@close_gui);

gui_data.curr_slice = 1;

% Layout (Using panel if available, or subplot fallback could be implemented, 
% but assuming panel based on Spine code preference)
gui_data.pp = panel();
gui_data.pp.pack('h', 2);
gui_data.pp.margin = [1 1 1 25]; 

% Controls string for Title
controls_str = ['\bfControls: \rm\leftarrow/\rightarrow: switch slice | ' ...
                'Click: add point | Backspace: delete last point | ' ...
                'C: clear all points | Space: toggle overlay | Enter: jump to slice | S: save'];
mode_str = upper(gui_data.mode);
gui_data.base_title = {['\bfControl point GUI (' mode_str '): \rmmatch points between sample and atlas'], controls_str};
gui_data.pp.title(gui_data.base_title);
gui_data.pp.fontname = 'Arial';

% Histology Axis
gui_data.histology_ax = gui_data.pp(1).select();
gui_data.histology_ax.YDir = 'reverse';
gui_data.histology_ax.Colormap = gray;
hold(gui_data.histology_ax, 'on'); axis(gui_data.histology_ax, 'image', 'off');

% Initial Histology Image
curr_image = volumeIdtoImage(gui_data.volume, gui_data.chooselist(gui_data.curr_slice, :));
if strcmp(mode, 'brain')
    % Brain mode specific: blank out background
    curr_image = blankImage_alt(curr_image, gui_data.chooselist(gui_data.curr_slice, 3:end));
end
curr_image = adapthisteq(curr_image);

gui_data.histology_im_h = imagesc(curr_image,...
    'Parent',gui_data.histology_ax,'ButtonDownFcn',@mouseclick_histology);
clim(gui_data.histology_ax, [0,255]);

% Overlay init
gui_data.histology_aligned_atlas_boundaries = ...
    plot(gui_data.histology_ax, nan, nan,...
    'r.','MarkerSize',3, 'PickableParts','none');

% Atlas Axis
gui_data.atlas_ax = gui_data.pp(2).select();
gui_data.atlas_ax.YDir = 'reverse';
gui_data.atlas_ax.Colormap = gray;
hold(gui_data.atlas_ax, 'on'); axis(gui_data.atlas_ax, 'image', 'off');

% Initial Atlas Slice
curr_atlas = volumeIdtoImage(gui_data.tv, gui_data.chooselist(gui_data.curr_slice, :));
gui_data.atlas_slice = gui_data.chooselist(gui_data.curr_slice, 1);

gui_data.atlas_im_h = imagesc(curr_atlas, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);
clim(gui_data.atlas_ax, [0,250]);

% Initialize Plot/Text placeholders
gui_data.h_pts_hist = plot(gui_data.histology_ax,nan,nan,'.g','MarkerSize',20);
gui_data.h_pts_atlas = plot(gui_data.atlas_ax,nan,nan,'.r','MarkerSize',20);
gui_data.h_text_hist = []; 
gui_data.h_text_atlas = []; 

% Initial Alignment logic
if isfield(gui_data,'histology_ccf_auto_alignment')
    gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;
end

guidata(gui_fig,gui_data);

% First Draw
align_ccf_to_histology(gui_fig);
update_slice(gui_fig); 

end

%==========================================================================
% Helper: Point Dimension Check
%==========================================================================
function points_cell = check_point_dims(points_cell)
    % Ensure points have 4 columns (x,y,z,timestamp). 
    % If 3 (legacy brain files), append NaN for time.
    if ~isempty(points_cell) && ~isempty(points_cell{1}) && size(points_cell{1}, 2) == 3
        points_cell = cellfun(@(x) [x, nan(size(x,1),1)], points_cell, 'UniformOutput', false);
    end
end

%==========================================================================
% GUI Callbacks (Unified)
%==========================================================================

function keypress(gui_fig,eventdata)
    gui_data = guidata(gui_fig);
    needs_update = false;

    switch eventdata.Key
        case 'return' % Go to slice
            input_slice = inputdlg(sprintf('Go to slice (max %d):',  size(gui_data.chooselist,1)));
            if ~isempty(input_slice)
                new_slice = str2double(input_slice{1});
                if ~isnan(new_slice) && new_slice >= 1 && new_slice <= size(gui_data.chooselist,1)
                    gui_data.curr_slice = new_slice;
                    needs_update = true;
                end
            end

        case 'leftarrow'
            gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
            needs_update = true;
            
        case 'rightarrow'
            gui_data.curr_slice = min(gui_data.curr_slice + 1,size(gui_data.chooselist,1));
            needs_update = true;
            
        case 'space' % Toggle overlay
            curr_vis = get(gui_data.histology_aligned_atlas_boundaries,'Visible');
            set(gui_data.histology_aligned_atlas_boundaries,'Visible', ...
                cell2mat(setdiff({'on','off'},curr_vis)));
            
        case 'c' % Clear current slice points
            gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,4);
            gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,4);
            needs_update = true;
            
        case 'backspace' % Delete last added point
            h_pts = gui_data.histology_control_points{gui_data.curr_slice};
            a_pts = gui_data.atlas_control_points{gui_data.curr_slice};
            
            % Find timestamps (use -inf if empty)
            t_h = -inf; t_a = -inf;
            if ~isempty(h_pts), t_h = h_pts(end, 4); end
            if ~isempty(a_pts), t_a = a_pts(end, 4); end
            
            if t_h > t_a
                gui_data.histology_control_points{gui_data.curr_slice}(end,:) = [];
            elseif t_a > -inf
                gui_data.atlas_control_points{gui_data.curr_slice}(end,:) = [];
            end
            
            needs_update = true;
            align_ccf_to_histology(gui_fig); 

        case 's' % Save
            atlas2histology_tform = gui_data.histology_ccf_manual_alignment;
            histology_control_points = gui_data.histology_control_points;
            atlas_control_points     = gui_data.atlas_control_points;
            
            save_fn = fullfile(gui_data.save_path, gui_data.save_filename);
            
            save(save_fn,'atlas2histology_tform', 'atlas_control_points', 'histology_control_points');
            
            % Feedback flash
            old_title = gui_data.base_title;
            gui_data.pp.title({'\bfSAVED!'});
            pause(0.2);
            gui_data.pp.title(old_title);
            disp(['Saved ' save_fn]);
    end

    if needs_update
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
    end
end

function mouseclick_histology(gui_fig,eventdata)
    gui_data = guidata(gui_fig);
    
    idim = gui_data.chooselist(gui_data.curr_slice, 2);
    alldims = 1:3;
    cpt = zeros(1,4);
    
    % Geometry
    cpt(alldims~=idim) = flip(eventdata.IntersectionPoint(1:2));
    cpt(idim) = gui_data.chooselist(gui_data.curr_slice, 1);
    
    cpt(4) = now; % Timestamp

    gui_data.histology_control_points{gui_data.curr_slice} = ...
        vertcat(gui_data.histology_control_points{gui_data.curr_slice}, cpt);

    guidata(gui_fig, gui_data);
    update_markers(gui_data, 'histology');
    check_and_align(gui_fig);
end

function mouseclick_atlas(gui_fig,eventdata)
    gui_data = guidata(gui_fig);

    idim = gui_data.chooselist(gui_data.curr_slice, 2);
    alldims = 1:3;
    cpt = zeros(1,4);
    
    % Geometry
    cpt(alldims~=idim) = flip(eventdata.IntersectionPoint(1:2));
    cpt(idim) = gui_data.atlas_slice;
    
    cpt(4) = now; % Timestamp

    gui_data.atlas_control_points{gui_data.curr_slice} = ...
        vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, cpt);

    guidata(gui_fig, gui_data);
    update_markers(gui_data, 'atlas');
    check_and_align(gui_fig);
end

function check_and_align(gui_fig)
    gui_data = guidata(gui_fig);
    nH = size(gui_data.histology_control_points{gui_data.curr_slice},1);
    nA = size(gui_data.atlas_control_points{gui_data.curr_slice},1);
    if nH == nA
        align_ccf_to_histology(gui_fig);
    end
end

%==========================================================================
% Plotting / Visuals
%==========================================================================

function update_markers(gui_data, type)
    % Efficiently redraw points and text numbers
    
    idim = gui_data.chooselist(gui_data.curr_slice, 2);
    alldims = 1:3;
    toplot = find(alldims~=idim);
    
    if strcmp(type, 'histology')
        ax = gui_data.histology_ax;
        h_plot = gui_data.h_pts_hist;
        pts = gui_data.histology_control_points{gui_data.curr_slice};
        old_text = gui_data.h_text_hist;
    else
        ax = gui_data.atlas_ax;
        h_plot = gui_data.h_pts_atlas;
        pts = gui_data.atlas_control_points{gui_data.curr_slice};
        old_text = gui_data.h_text_atlas;
    end
    
    % 1. Update Scatter Plot
    if ~isempty(pts)
        set(h_plot, 'XData', pts(:,toplot(2)), 'YData', pts(:,toplot(1)));
    else
        set(h_plot, 'XData', nan, 'YData', nan);
    end
   
    % 2. Update Text Numbers
    if isa(old_text, 'matlab.graphics.Graphics') 
         delete(old_text(isvalid(old_text))); 
    elseif ~isempty(old_text) && isnumeric(old_text)
         % Fallback
    end

    new_text = gobjects(size(pts,1), 1);
    if ~isempty(pts)
        x_coord = pts(:,toplot(2));
        y_coord = pts(:,toplot(1));
        for i = 1:size(pts,1)
            new_text(i) = text(ax, x_coord(i), y_coord(i), num2str(i), ...
                'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
                'PickableParts', 'none');
        end
    end
    
    if strcmp(type, 'histology')
        gui_data.h_text_hist = new_text;
    else
        gui_data.h_text_atlas = new_text;
    end
    guidata(gui_data.pp.figure, gui_data);
end

function update_slice(gui_fig)
    gui_data = guidata(gui_fig);

    tform  = affinetform3d(gui_data.histology_ccf_manual_alignment);
    idim   = gui_data.chooselist(gui_data.curr_slice, 2);
    induse = gui_data.chooselist(gui_data.curr_slice, 1);
    
    % Get Histology Image
    curr_image = volumeIdtoImage(gui_data.volume, gui_data.chooselist(gui_data.curr_slice, :));
    curr_image = adapthisteq(curr_image);
    if strcmp(gui_data.mode, 'brain')
         curr_image = blankImage_alt(curr_image, gui_data.chooselist(gui_data.curr_slice, 3:end));
    end
    
    % Auto-calculate matching atlas slice
    [~,iy,ix] = blankImage_alt(curr_image, gui_data.chooselist(gui_data.curr_slice, 3:end));
    [xx,yy]   = meshgrid(1:size(curr_image,2), 1:size(curr_image,1));
    yy = yy(iy, ix); xx = xx(iy, ix);
    itform = tform.invert;
    
    switch idim
        case 1
            [~,~,zn] = tform.transformPointsInverse(yy(:),induse*ones(numel(xx),1),xx(:));
            [~, yl, ~] = itform.outputLimits([1 size(curr_image,1)], [induse induse], [1 size(curr_image,2)]);
            sluse = round(median(yl));
        case 2
            [~,~,zn] = tform.transformPointsInverse(induse*ones(numel(xx),1),yy(:),xx(:));
            [xl, ~, ~] = itform.outputLimits([induse induse], [1 size(curr_image,1)], [1 size(curr_image,2)]);
            sluse = round(median(xl));
        case 3
            [~,~,zn] = tform.transformPointsInverse(xx(:),yy(:),induse*ones(numel(xx),1));
            [xl, yl, zl] = itform.outputLimits([1 size(curr_image,2)], [1 size(curr_image,1)], [induse induse]);
            sluse = round(median(zl));
    end
    sluse = max(sluse, 1);

    % Use manual points context if available
    cpointsatlas = gui_data.atlas_control_points{gui_data.curr_slice};
    if ~isempty(cpointsatlas)
        gui_data.atlas_slice = round(median(cpointsatlas(:,idim)));
    else
        gui_data.atlas_slice = sluse;
    end

    % Update Histology Plot
    currlim = getImageLimits(curr_image, 0.001);
    set(gui_data.histology_im_h,'CData', curr_image);
    gui_data.histology_ax.CLim = [0 currlim(2)];

    title(gui_data.histology_ax, sprintf('Sample Slice %d', induse), 'FontSize', 12);
    set(gui_data.histology_aligned_atlas_boundaries, 'XData',nan, 'YData',nan);

    guidata(gui_fig, gui_data);
    update_markers(gui_data, 'histology');
    align_ccf_to_histology(gui_fig);
    update_atlas_slice(gui_fig);
end

function update_atlas_slice(gui_fig)
    gui_data = guidata(gui_fig);

    idim = gui_data.chooselist(gui_data.curr_slice, 2);
    sluse = gui_data.atlas_slice;
    atlasize = size(gui_data.tv);
    
    sluse = max(min(sluse, atlasize(idim)), 1); 
    
    curr_atlas = volumeIdtoImage(gui_data.tv, [sluse idim]);
    curr_atlas = adapthisteq(curr_atlas);

    set(gui_data.atlas_im_h,'CData', curr_atlas);
    gui_data.atlas_ax.Title.String = sprintf("Atlas slice %d/%d", sluse, atlasize(idim));

    update_markers(gui_data, 'atlas');
end

function scroll_atlas_slice(gui_fig,eventdata)
    gui_data = guidata(gui_fig);
    gui_data.atlas_slice = gui_data.atlas_slice + 2* eventdata.VerticalScrollCount;
    guidata(gui_fig, gui_data);
    update_atlas_slice(gui_fig);
end

%==========================================================================
% Alignment Logic (Unified)
%==========================================================================

function align_ccf_to_histology(gui_fig)
    gui_data = guidata(gui_fig);

    Nmin = 16;
    % Extract points (cols 1-3, ignore timestamp)
    cptsatlas = cat(1, gui_data.atlas_control_points{:});
    if ~isempty(cptsatlas), cptsatlas = cptsatlas(:, 1:3); end
    
    cptshistology = cat(1, gui_data.histology_control_points{:});
    if ~isempty(cptshistology), cptshistology = cptshistology(:, 1:3); end
    
    % Rearrange dimensions [2 1 3] for transform
    if ~isempty(cptsatlas), cptsatlas = cptsatlas(:, [2 1 3]); end
    if ~isempty(cptshistology), cptshistology = cptshistology(:, [2 1 3]); end

    tform = affinetform3d; 

    if size(cptshistology,1) == size(cptsatlas,1) && ...
            (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
        
        [tform, mse] = fitAffineTrans3D(cptsatlas, cptshistology);

        status_str = sprintf('\\bfCurrent Fit: \\rmMSE: %2.2f | Npoints: %d', mse, size(cptshistology,1));
        gui_data.pp.title([gui_data.base_title, {status_str}]);
        
        gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rvolume);
        
    elseif isfield(gui_data,'histology_ccf_auto_alignment')
        tform = affinetform3d(gui_data.histology_ccf_auto_alignment);
        gui_data.pp.title([gui_data.base_title, {'\bfCurrent Fit: \rmUsing initial auto-alignment'}]);
        gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rvolume);
    else
        gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rvolume);
        gui_data.pp.title([gui_data.base_title, {'\bfCurrent Fit: \rmIdentity/No alignment'}]);
    end

    % Draw boundaries
    curr_slice_warp = volumeIdtoImage(gui_data.volwrap, gui_data.chooselist(gui_data.curr_slice,:));
    av_warp_boundaries = round(conv2(curr_slice_warp,ones(3)./9,'same')) ~= curr_slice_warp;
    [row,col] = ind2sub(size(curr_slice_warp), find(av_warp_boundaries));

    set(gui_data.histology_aligned_atlas_boundaries, 'XData', col, 'YData', row);

    gui_data.histology_ccf_manual_alignment = tform.A;
    guidata(gui_fig, gui_data);
end

function close_gui(gui_fig,~)
    gui_data = guidata(gui_fig);
    user_confirm = questdlg('Save changes?','Confirm exit', 'Yes', 'No', 'Cancel', 'Yes');
    switch user_confirm
        case 'Yes'
            atlas2histology_tform = gui_data.histology_ccf_manual_alignment;
            histology_control_points = gui_data.histology_control_points;
            atlas_control_points = gui_data.atlas_control_points;
            
            save_fn = fullfile(gui_data.save_path, gui_data.save_filename);
            save(save_fn,'atlas2histology_tform', 'atlas_control_points', 'histology_control_points');
            disp(['Saved ' save_fn]);
            delete(gui_fig);
        case 'No'
            delete(gui_fig);
    end   
end