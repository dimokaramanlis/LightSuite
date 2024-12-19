function reorder_dimensions(~,~,histology_toolbar_gui)

% Get images (from path in GUI)
histology_toolbar_guidata = guidata(histology_toolbar_gui);


volume_dir = dir(fullfile(histology_toolbar_guidata.save_path,'*.tif'));
volpath    = fullfile(volume_dir.folder, volume_dir.name);
voldown    = readDownStack(volpath);


% Plot all (downsampled) images
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.2; % width/length
gui_width_fraction = 0.6; % fraction of screen width to occupy
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... % left x
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... % bottom y
    gui_width_px,gui_width_px/gui_aspect_ratio]; % width, height

gui_fig = figure('Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position);
tile_h = tiledlayout('vertical','TileSpacing',"tight");
image_h = gobjects(3,1);
for curr_slice = 1:3
    nexttile; 
    
    Nall    = size(voldown, curr_slice);
    idsplot = round(linspace(Nall*0.2, Nall*0.8, 5));
    S.type = '()';
    S.subs = repmat({':'}, 1, 3);  % create a cell array with ':' for all dimensions
    S.subs{curr_slice} = idsplot;             % replace the id-th dimension with the index_value
    result = subsref(voldown, S);
    avdims = 1:3;
    avdims(curr_slice) = [];
    image_h(curr_slice) = montage(permute(result, [avdims curr_slice]),'Size',[1 5]);
    
    % imagesc(imresize(slice_im{curr_slice},1/10,'nearest'));
    axis image off;
end

% Set click function
[image_h.ButtonDownFcn] = deal({@click_slice,gui_fig});

% Title with directions
title(tile_h,{'Click to assign/un-assign dimension order', ...
    'How to set: 1-Coronal, 2-Axial, 3-Saggital'},...
    'FontSize',12);

% Package image handles and slice number index in figure
gui_data = struct;
gui_data.slice_idx = nan(3, 1);
gui_data.image_h   = image_h;
gui_data.volume    = voldown;
gui_data.volpath   = volpath;
guidata(gui_fig,gui_data)

end

function click_slice(obj,eventdata,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Get current index of slice
curr_ax = find(gui_data.image_h == obj);

if ~any(gui_data.slice_idx == curr_ax)
    % If slice isn't assigned, assigned next number

    % Get number to assign
    curr_idx_assign = find(isnan(gui_data.slice_idx),1);

    % Assign number to currently ordered slice
    gui_data.slice_idx(curr_idx_assign) = curr_ax;

    % Write number on axis
    text(get(obj,'parent'),20,20,num2str(curr_idx_assign),'FontSize',20,'Color','w')

else
    % If slice is already assigned, remove assignment
    gui_data.slice_idx(gui_data.slice_idx == curr_ax) = NaN;

    % Clear number from axis 
    delete(findobj(get(obj,'parent'),'type','text'));
end

% Upload gui data
guidata(gui_fig,gui_data)

% If all slices assigned, close and save
if ~any(isnan(gui_data.slice_idx))
    close(gui_fig);
    save_reordered_volume(gui_data);
end

end

function save_reordered_volume(gui_data)

% Make re-ordered filenames (with '_reorder to avoid overwriting)
volsave = permute(gui_data.volume, gui_data.slice_idx);
reorder_target = strrep(gui_data.volpath,'.tif','_reorder.tif');
saveLightsheetVolume(volsave, reorder_target)
movefile(reorder_target, gui_data.volpath);
disp('Saved re-ordered volume');

end





