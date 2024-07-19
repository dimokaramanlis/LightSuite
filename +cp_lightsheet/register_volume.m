function register_volume(~,~,histology_toolbar_gui, dimalign)
% Part of cp_lightsheet toolbox
%
% Choose CCF atlas slices corresponding to histology slices

% Initialize guidata
gui_data = struct;
% Store toolbar handle
gui_data.histology_toolbar_gui = histology_toolbar_gui;

% Load atlas
% allen_atlas_path = fileparts(which('template_volume_10um.npy'));
% if isempty(allen_atlas_path)
%     error('No CCF atlas found (add CCF atlas to path)')
% end
% disp('Loading Allen CCF atlas...')
% 
% gui_data.tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
% gui_data.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
% gui_data.st = cp_lightsheet.loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

[tv, av, atlas_path] = readGubraAtlas();

gui_data.tv = tv;
gui_data.av = av;
gui_data.st = cp_lightsheet.loadStructureTree(fullfile(atlas_path,'structure_tree_safe_2017.csv'));
disp('Done.')

% Get images (from path in GUI)
histology_toolbar_guidata = guidata(histology_toolbar_gui);

volume_dir     = dir(fullfile(histology_toolbar_guidata.save_path,'*.tif'));
volpath        = fullfile(volume_dir.folder, volume_dir.name);
gui_data.volume = readDownStack(volpath);

[Ny, Nx, Nz ] = size(gui_data.volume);
atlsize       = size(gui_data.av);

S.type = '()';
S.subs = repmat({':'}, 1, 3);  % create a cell array with ':' for all dimensions
switch dimalign
    case 'coronal'
        curr_slice = 1;
    case 'axial'
        curr_slice = 2;
    case 'sagittal'
        curr_slice = 3;
end
S.subs{curr_slice} = 1;             % replace the id-th dimension with the index_value
gui_data.dim_index = curr_slice;
gui_data.curr_histology_slice = 1;
gui_data.atlas_slice_point = [1 1 1];
ikeep = 1:3;
ikeep = ~(ikeep == curr_slice);

% Set save path (from toolbar GUI
gui_data.save_path = histology_toolbar_guidata.save_path;

% Create figure, set button functions
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7; % width/length
gui_width_fraction = 0.6; % fraction of screen width to occupy
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... % left x
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... % bottom y
    gui_width_px,gui_width_px/gui_aspect_ratio]; % width, height

gui_fig = figure('WindowScrollWheelFcn',@scroll_atlas_slice, ...
    'KeyPressFcn',@keypress,'Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position, ...
    'CloseRequestFcn',@close_gui);

% p = panel(gui_fig,'add');
% p.pack('h', 2);
% p.de.margin = 10;
% p.margin = [1 15 15 15];
% Set up axis for histology image
gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
gui_data.histology_ax.Colormap = gray;
hold on; axis image off;
gui_data.histology_im_h = image(squeeze(subsref(gui_data.volume , S)),'Parent',gui_data.histology_ax);
title(gui_data.histology_ax,'No saved atlas position');

atlaskeep = atlsize(ikeep);
% Set up 3D atlas axis
gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse', ...
    'XTick',[1, atlaskeep(2)],'XTickLabel',{'Front','Back'}, ...
    'YTick',[1, atlaskeep(1)],'YTickLabel',{'Left','Right'});
hold on
axis  equal manual
xlim([1,atlaskeep(2)]);
ylim([1,atlaskeep(1)]);
colormap(gui_data.atlas_ax,'gray');
clim([0,400]);
gui_data.atlas_title = title(sprintf('Slice position: %d',0));

% Create CCF colormap
% (copied from cortex-lab/allenCCF/setup_utils
ccf_color_hex = gui_data.st.color_hex_triplet;
ccf_color_hex(cellfun(@numel,ccf_color_hex)==5) = {'019399'}; % special case where leading zero was evidently dropped
ccf_cmap_c1 = cellfun(@(x)hex2dec(x(1:2)), ccf_color_hex, 'uni', false);
ccf_cmap_c2 = cellfun(@(x)hex2dec(x(3:4)), ccf_color_hex, 'uni', false);
ccf_cmap_c3 = cellfun(@(x)hex2dec(x(5:6)), ccf_color_hex, 'uni', false);
gui_data.ccf_cmap = ...
    horzcat(vertcat(ccf_cmap_c1{:}),vertcat(ccf_cmap_c2{:}),vertcat(ccf_cmap_c3{:}))./255;

% Set mode for atlas view (can be either TV, AV, or TV-AV)
gui_data.atlas_mode = 'TV';

% Create slice object and first slice point
imstart = volumeIdtoImage(gui_data.tv, [gui_data.atlas_slice_point gui_data.dim_index]);
gui_data.atlas_slice_plot = imagesc(gui_data.atlas_ax, imstart); % Slice on 3D atlas

% Set up atlas parameters to save for histology
gui_data.slice_points = nan(size(gui_data.volume, curr_slice),1);

% Upload gui data
guidata(gui_fig,gui_data);

% Draw the first slice
update_atlas_slice(gui_fig);

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Controls: \rm' ...
    'Left/right arrows: cycle through volume dimension', ...
    'Shift + arrows: change atlas rotation', ...
    'm : change atlas display mode (TV/AV/TV-AV overlay)', ...
    'Scroll wheel: move CCF slice in/out of plane', ...
    'Enter: set current histology and CCF slice pair'}, ...
    'Controls',CreateStruct);

end 

function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

shift_on = any(strcmp(eventdata.Modifier,'shift'));


dim_index = gui_data.dim_index;

switch eventdata.Key
    
    % Left/right: cycle through histology slices
    % (if there's a saved plane point, move atlas to that position)
    % Shift + arrow keys: rotate atlas slice
    case 'leftarrow'
        if ~shift_on
            new_slice_index = max(gui_data.curr_histology_slice - 1, 1);
            gui_data.curr_histology_slice = new_slice_index;
            guidata(gui_fig,gui_data);
            update_histology_slice(gui_fig);
        elseif shift_on
            set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [1,0]);
            update_atlas_slice(gui_fig)
        end
    case 'rightarrow'
        if ~shift_on
            max_index = size(gui_data.volume, dim_index);
            new_slice_index = min(gui_data.curr_histology_slice + 1, max_index);
            gui_data.curr_histology_slice = new_slice_index;
            guidata(gui_fig,gui_data);
            update_histology_slice(gui_fig);
        elseif shift_on
            set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
            update_atlas_slice(gui_fig)
        end
    case 'uparrow'
        if shift_on
            set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,-1]);
            update_atlas_slice(gui_fig)
        end
    case 'downarrow'
        if shift_on
            set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,1]);
            update_atlas_slice(gui_fig)
        end

    % M key: switch atlas display mode
    case 'm'
        atlas_slice_modes = {'TV','AV','TV-AV'};
        curr_atlas_mode_idx = strcmp(gui_data.atlas_mode,atlas_slice_modes);
        gui_data.atlas_mode = atlas_slice_modes{circshift(curr_atlas_mode_idx,[0,1])};
        guidata(gui_fig,gui_data);
        update_atlas_slice(gui_fig);

    % Enter: save slice coordinates
    case 'return'        
        % Store camera vector and point
        % (Note: only one camera vector used for all slices, overwrites)
        gui_data.slice_points(gui_data.curr_histology_slice) = ...
            gui_data.atlas_slice_point(dim_index);
        guidata(gui_fig,gui_data);
                
        update_histology_slice(gui_fig);
        title(gui_data.histology_ax,'New saved atlas position');
        
end

end

function update_histology_slice(gui_fig)
% Draw histology slice (and move atlas if saved position)

% Get guidata
gui_data = guidata(gui_fig);
dim_index = gui_data.dim_index;


S.type = '()';
S.subs = repmat({':'}, 1, 3);  % create a cell array with ':' for all dimensions
S.subs{dim_index} = gui_data.curr_histology_slice;             % replace the id-th dimension with the index_value


displayed_slice = squeeze(subsref(gui_data.volume, S));

% Set next histology slice
set(gui_data.histology_im_h,'CData',displayed_slice)

% If there's a saved atlas position, move atlas to there
if all(~isnan(gui_data.slice_points(gui_data.curr_histology_slice,:)))
    gui_data.atlas_slice_point(dim_index) = ...
        gui_data.slice_points(gui_data.curr_histology_slice,:);
    title(gui_data.histology_ax,'Saved atlas position')
    guidata(gui_fig,gui_data);
    update_atlas_slice(gui_fig);
else
    title(gui_data.histology_ax,'No saved atlas position')
end

% Upload gui data
guidata(gui_fig, gui_data);

end

function cam_vector = get_camera_vector(gui_data)
% Get the camera viewing vector to define atlas slice plane

% Grab current camera angle
% (normalized line from the camera to the center)
curr_campos = campos(gui_data.atlas_ax);
curr_camtarget = camtarget(gui_data.atlas_ax);
cam_vector = (curr_camtarget - curr_campos)./norm(curr_camtarget - curr_campos);

end

function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point
gui_data.atlas_slice_point(gui_data.dim_index) = ...
    gui_data.atlas_slice_point(gui_data.dim_index) + eventdata.VerticalScrollCount*2;

gui_data.atlas_slice_point = max(gui_data.atlas_slice_point, 1);
gui_data.atlas_slice_point = min(gui_data.atlas_slice_point, size(gui_data.tv));

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_atlas_slice(gui_fig)

end



function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

% Get slice (larger spacing for faster pulling)
[tv_slice,av_slice,slice_coords] = grab_atlas_slice(gui_data);

% Update the slice display (depending on display mode)
switch gui_data.atlas_mode
    case 'TV'
        atlas_slice = tv_slice;
        colormap(gui_data.atlas_ax, gray);clim([0,516]);
    case 'AV'
        av_boundaries = round(conv2(av_slice,ones(2)./4,'same')) ~= av_slice;
        atlas_slice = imoverlay(mat2gray(tv_slice,[0,516]),av_boundaries,'r');
        clim([0,1]);
    case 'TV-AV'
        atlas_slice = av_slice;
        colormap(gui_data.atlas_ax, gui_data.ccf_cmap)
        clim(gui_data.atlas_ax,[1,size(gui_data.ccf_cmap,1)])
end
set(gui_data.atlas_slice_plot,'CData',atlas_slice);
gui_data.atlas_ax.XLim = [0.5 size(tv_slice,2)];
gui_data.atlas_ax.YLim = [0.5 size(tv_slice,1)];

% Upload gui_data
guidata(gui_fig, gui_data);

end

function [tv_slice,av_slice,slice_coords] = grab_atlas_slice(gui_data)
% Grab anatomical and labelled atlas within slice

idim = gui_data.dim_index;
tv_slice = volumeIdtoImage(gui_data.tv, [gui_data.atlas_slice_point(idim) idim]);
av_slice = volumeIdtoImage(gui_data.av, [gui_data.atlas_slice_point(idim) idim]);

slice_coords = size(gui_data.volume)/2;
slice_coords(gui_data.dim_index) = gui_data.atlas_slice_point(idim);


% Update slice position title
plane_offset_mm = gui_data.atlas_slice_point(idim)/100; % CCF = 10um voxels
set(gui_data.atlas_title,'string', ...
    sprintf('Slice position: %.2f mm',plane_offset_mm));

end

function close_gui(gui_fig,~)

% Get guidata
gui_data = guidata(gui_fig);

opts.Default = 'Yes';
opts.Interpreter = 'tex';
user_confirm = questdlg('\fontsize{14} Save?','Confirm exit',opts);
switch user_confirm
    case 'Yes'

        indsassign = find(any(~isnan(gui_data.slice_points),2));
        save_fn = fullfile(gui_data.save_path,'histology_ccf.mat');
        indids        = [];
        histology_ccf = [];
        if exist(save_fn, 'file')
            user_confirm2 = questdlg('\fontsize{14} Keep old slices?','Confirm exit',opts);
            switch user_confirm2
                case 'Yes'
                    oldslices     = load(save_fn);
                    indids        = oldslices.indids;
                    histology_ccf = oldslices.histology_ccf;
            end
        end

        % Go through each slice, pull full-resolution atlas slice and
        % corrsponding coordinates
        
        histology_ccf_init = cell(size(indsassign, 1),1);
        histology_ccf_curr = struct( ...
            'tv_slices',histology_ccf_init, ...
            'slice_coords',histology_ccf_init);

        h = waitbar(0,'Saving atlas slices...');
        for curr_slice = 1:numel(indsassign)
            gui_data.atlas_slice_point(gui_data.dim_index) = gui_data.slice_points( indsassign(curr_slice),:);
            [histology_ccf_curr(curr_slice).tv_slices, ...
                ~, histology_ccf_curr(curr_slice).slice_coords] = ...
                grab_atlas_slice(gui_data);
            waitbar(curr_slice/numel(indsassign),h, ...
                ['Saving atlas slices (' num2str(curr_slice) '/' num2str(numel(indsassign)) ')...']);
        end
        close(h);

        indidsnew     = cat(2, indsassign, gui_data.dim_index*ones(size(indsassign)));
        indids        = cat(1, indids, indidsnew);
        histology_ccf = cat(1, histology_ccf, histology_ccf_curr);
        save(save_fn,'histology_ccf','indids','-v7.3');
        delete(gui_fig);

    case 'No'
        % Close without saving
        delete(gui_fig);

    case 'Cancel'
        % Do nothing

end 

% Update toolbar GUI
cp_lightsheet.update_toolbar_gui(gui_data.histology_toolbar_gui);

end











