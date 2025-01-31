function matchControlPoints(opts)
% Part of AP_histology toolbox
%
% Manually align histology slices and matched CCF slices

% Initialize guidata
gui_data = struct;

% Load atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
factv       = 255/single(max(gui_data.tv,[],"all"));
gui_data.tv = uint8(single(gui_data.tv)*factv);
gui_data.tv = imresize3(gui_data.tv,opts.downfac_reg);
gui_data.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
gui_data.av = imresize3(gui_data.av,opts.downfac_reg, "Method","nearest");
gui_data.st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
gui_data.Rmoving  = imref3d(size(gui_data.av));
disp('Done.')

gui_data.save_path = opts.savepath;

volume_dir      = dir(fullfile(opts.savepath,'sample_register*.tif'));
volpath         = fullfile(volume_dir.folder, volume_dir.name);
volload         = readDownStack(volpath);
volload         = permute(volload, opts.permute_sample_to_atlas);
Rvolume         = imref3d(size(volload));
gui_data.volume = imwarp(volload, Rvolume, opts.original_trans, 'OutputView',gui_data.Rmoving);
factv       = 255/single(max(gui_data.volume,[],"all"));
gui_data.volume = uint8(single(gui_data.volume)*factv);
chooselist      = generate_cp_list(gui_data.volume);
gui_data.chooselist = chooselist;

% chooselist = cell(3,1);
% for ii = 1:3
%     sids = round(size(gui_data.volume, ii)*0.1):round(size(gui_data.volume, ii)*0.9);
%     chooselist{ii} = [sids' ii * ones(numel(sids),1) ];
% end
% chooselist = cat(1, chooselist{:});
% rng(1);
% gui_data.chooselist = chooselist(randperm(size(chooselist, 1)),:);


% Load corresponding CCF slices
% ccf_slice_fn = fullfile(gui_data.save_path,'histology_ccf.mat');
% aldata       = load(ccf_slice_fn);

% Nslices = size(aldata.indids, 1);
% rng(1);
% irand   = randperm(Nslices);
% gui_data.histology_ccf = aldata.histology_ccf(irand);
% gui_data.volindids     = aldata.indids(irand,:);



% Load automated alignment
auto_ccf_alignment_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
if exist(auto_ccf_alignment_fn,'file')
    oldtform = load(auto_ccf_alignment_fn);
    gui_data.histology_ccf_auto_alignment = oldtform.atlas2histology_tform;
end

% Create figure, set button functions
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7; % width/length
gui_width_fraction = 0.6; % fraction of screen width to occupy
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... % left x
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... % bottom y
    gui_width_px,gui_width_px/gui_aspect_ratio]; % width, height

gui_fig = figure('KeyPressFcn',@keypress, ...
    'WindowScrollWheelFcn',@scroll_atlas_slice,...
    'Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position, ...
    'CloseRequestFcn',@close_gui);


% gui_data.curr_slice = randperm(numel(chooselist), 1);
% curr_image = volumeIdtoImage(gui_data.volume, gui_data.volindids(1, :));
gui_data.curr_slice = 1;
curr_image = volumeIdtoImage(gui_data.volume, chooselist(gui_data.curr_slice, :));
curr_image = blankImage(curr_image, chooselist(gui_data.curr_slice, 3:end));
curr_image = adapthisteq(curr_image);

% Set up axis for histology image
gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
set(gui_data.histology_ax,'Position',[0,0,0.5,0.9]);
hold on; colormap(gray); axis image off;
gui_data.histology_im_h = imagesc(curr_image,...
    'Parent',gui_data.histology_ax,'ButtonDownFcn',@mouseclick_histology);
clim([0,200]);
% Set up histology-aligned atlas overlay
% (and make it invisible to mouse clicks)
% histology_aligned_atlas_boundaries_init = ...
%     zeros(size(curr_image));
% gui_data.histology_aligned_atlas_boundaries = ...
%     imagesc(histology_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
%     'AlphaData',histology_aligned_atlas_boundaries_init,'PickableParts','none');

histology_aligned_atlas_boundaries_init = ...
    zeros(size(curr_image));
gui_data.histology_aligned_atlas_boundaries = ...
    plot(histology_aligned_atlas_boundaries_init(:,1), histology_aligned_atlas_boundaries_init(:,2),...
    'r.','MarkerSize',3, 'Parent',gui_data.histology_ax, 'PickableParts','none');

curr_atlas = volumeIdtoImage(gui_data.tv, chooselist(gui_data.curr_slice, :));
gui_data.atlas_slice = chooselist(gui_data.curr_slice, 1);

% Set up axis for atlas slice
gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse'); 
set(gui_data.atlas_ax,'Position',[0.5,0,0.5,0.9]);
hold on; axis image off; colormap(gray); clim([0,250]);
gui_data.atlas_im_h = imagesc(curr_atlas, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);

% Initialize alignment control points and tform matricies
gui_data.histology_control_points = repmat({zeros(0,3)},size(gui_data.chooselist, 1),1);
gui_data.atlas_control_points     = repmat({zeros(0,3)},size(gui_data.chooselist, 1),1);

gui_data.histology_control_points_plot = plot(gui_data.histology_ax,nan,nan,'.g','MarkerSize',20);
gui_data.atlas_control_points_plot = plot(gui_data.atlas_ax,nan,nan,'.r','MarkerSize',20);

% If there was previously auto-alignment, intitialize with that
if isfield(gui_data,'histology_ccf_auto_alignment')
    gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;
end

% Upload gui data
guidata(gui_fig,gui_data);

% Initialize alignment
align_ccf_to_histology(gui_fig);

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Controls: \rm' ...
    'Left/right: switch slice' ...
    'click: set reference points for manual alignment (3 minimum)', ...
    'space: toggle alignment overlay visibility', ...
    'c: clear reference points', ...
    's: save'}, ...
    'Controls',CreateStruct);

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % left/right arrows: move slice
    case 'leftarrow'
        gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    case 'rightarrow'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice + 1,size(gui_data.chooselist,1));
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % space: toggle overlay visibility
    case 'space'
        curr_visibility = ...
            get(gui_data.histology_aligned_atlas_boundaries,'Visible');
        set(gui_data.histology_aligned_atlas_boundaries,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
        
    % c: clear current points
    case 'c'
        gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,3);
        gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,3);
        
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % s: save
    case 's'
        atlas2histology_tform = ...
            gui_data.histology_ccf_manual_alignment{end};
        histology_control_points = gui_data.histology_control_points;
        atlas_control_points     = gui_data.atlas_control_points;
        save_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
        save(save_fn,'atlas2histology_tform', 'atlas_control_points', 'histology_control_points');
        disp(['Saved ' save_fn]);
        
end

end


function mouseclick_histology(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);


idim = gui_data.chooselist(gui_data.curr_slice, 2);
alldims = 1:3;
cpt     = zeros(1,3);
cpt(alldims~=idim) =  flip(eventdata.IntersectionPoint(1:2));
cpt(idim) = gui_data.chooselist(gui_data.curr_slice, 1);
toplot = find(alldims~=idim);


% Add clicked location to control points
gui_data.histology_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.histology_control_points{gui_data.curr_slice}, ...
    cpt);

set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(1)));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 5 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 5)
    align_ccf_to_histology(gui_fig)
end

end


function mouseclick_atlas(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

idim = gui_data.chooselist(gui_data.curr_slice, 2);
dval = gui_data.atlas_slice;
alldims = 1:3;
cpt     = zeros(1,3);
cpt(alldims~=idim) =  flip(eventdata.IntersectionPoint(1:2));
cpt(idim) = dval;

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    cpt);

toplot = find(alldims~=idim);
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 5 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 5)
    align_ccf_to_histology(gui_fig)
end

end


function align_ccf_to_histology(gui_fig)

% Get guidata 
gui_data = guidata(gui_fig);


Nmin = 8;
cptsatlas     = cat(1, gui_data.atlas_control_points{:});
cptsatlas     = cptsatlas(:, [2 1 3]);
cptshistology = cat(1, gui_data.histology_control_points{:});
cptshistology = cptshistology(:, [2 1 3]);

if size(cptshistology,1) == size(cptsatlas,1) && ...
        (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
    
    % If same number of >= 3 control points, use control point alignment
    [tform, mse] = fitAffineTrans3D(cptsatlas, cptshistology);
    % tform = estgeotform3d(cptsatlas, cptshistology, 'similarity');
    gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rmoving);

    title(gui_data.histology_ax, sprintf('New alignment, mse = %2.2f, Npoints = %d', ...
        mse, size(cptshistology,1)));


elseif size(gui_data.histology_control_points{gui_data.curr_slice},1) >= 1 ||  ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) >= 1
    % If less than 3 or nonmatching points, use auto but don't draw
    title(gui_data.histology_ax,'New alignment');

    % Upload gui data
    guidata(gui_fig, gui_data);
    return

elseif isfield(gui_data,'histology_ccf_auto_alignment')
    % If no points, use automated outline if available
    tform = gui_data.histology_ccf_auto_alignment;
    tform = affinetform3d(tform);
    title(gui_data.histology_ax,'Previous alignment');
    gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rmoving);
else
    % If nothing available, use identity transform
    % allpts = ones(size(gui_data.volindids,1), 3);
    % allpts = allpts.*size(gui_data.volume)/2;
    % 
    % allatlaspts = cat(1,gui_data.histology_ccf(:).slice_coords);
    % 
    % for idim = 1:3
    %     icurr = gui_data.volindids(:,2) == idim;
    %     allpts(icurr, idim) = gui_data.volindids(icurr,1);
    % end
    tform = affinetform3d;
    % allatlaspts = allatlaspts(:, [2 1 3]);
    % allpts = allpts(:, [2 1 3]);
    % tform = fitAffineTrans3D(allatlaspts, allpts);

    gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'OutputView',gui_data.Rmoving);
end


curr_slice_warp    = volumeIdtoImage(gui_data.volwrap, gui_data.chooselist(gui_data.curr_slice,:));
av_warp_boundaries = round(conv2(curr_slice_warp,ones(3)./9,'same')) ~= curr_slice_warp;

[row,col] = ind2sub(size(curr_slice_warp), find(av_warp_boundaries));


 % set(gui_data.histology_aligned_atlas_boundaries, ...
 %    'CData',av_warp_boundaries, ...
 %    'AlphaData',av_warp_boundaries*0.6);
set(gui_data.histology_aligned_atlas_boundaries, ...
'XData', col, 'YData', row);


% Update transform matrix
gui_data.histology_ccf_manual_alignment = tform.A;

% Upload gui data
guidata(gui_fig, gui_data);

end


function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

tform = affinetform3d(gui_data.histology_ccf_manual_alignment);

idim    = gui_data.chooselist(gui_data.curr_slice, 2);
alldims = 1:3;
toplot  = find(alldims~=idim);

% Set next histology slice
curr_image = volumeIdtoImage(gui_data.volume, gui_data.chooselist(gui_data.curr_slice, :));
[curr_image,iy,ix] = blankImage(curr_image,  gui_data.chooselist(gui_data.curr_slice, 3:end));

induse     = gui_data.chooselist(gui_data.curr_slice, 1);
[xx,yy]    = meshgrid(1:size(curr_image,2), 1:size(curr_image,1));
itform     = tform.invert;
yy = yy(iy, ix);
xx = xx(iy, ix);
switch idim
    case 1
        [xn,yn,zn]   = tform.transformPointsInverse(yy(:),induse*ones(numel(xx),1),xx(:));
        [xl, yl, zl] = itform.outputLimits([1 size(curr_image,1)], [induse induse], [1 size(curr_image,2)]);
        sluse = round(median(yl));
    case 2
        [xn,yn,zn] = tform.transformPointsInverse(induse*ones(numel(xx),1),yy(:),xx(:));
        [xl, yl, zl] = itform.outputLimits([induse induse], [1 size(curr_image,1)], [1 size(curr_image,2)]);
        sluse = round(median(xl));
    case 3
        [xn,yn,zn] = tform.transformPointsInverse(xx(:),yy(:),induse*ones(numel(xx),1));
        [xl, yl, zl] = itform.outputLimits([1 size(curr_image,2)], [1 size(curr_image,1)], [induse induse]);
        sluse = round(median(zl));
end

sluse = max(sluse, 1);
gui_data.atlas_slice = sluse;


curr_image = adapthisteq(curr_image);
set(gui_data.histology_im_h,'CData',curr_image)
% Plot control points for slice
set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(1)));

histology_aligned_atlas_boundaries_init = nan(1,2);
set(gui_data.histology_aligned_atlas_boundaries, ...
    'XData',histology_aligned_atlas_boundaries_init(:,1), 'YData',histology_aligned_atlas_boundaries_init(:,2));



% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

update_atlas_slice(gui_fig)



end


function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point
gui_data.atlas_slice = ...
    gui_data.atlas_slice + eventdata.VerticalScrollCount;

idim = gui_data.chooselist(gui_data.curr_slice, 2);

gui_data.atlas_slice = max(gui_data.atlas_slice, 1);
gui_data.atlas_slice = min(gui_data.atlas_slice, size(gui_data.tv, idim));

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_atlas_slice(gui_fig)

end

function close_gui(gui_fig,~)

% Get guidata
gui_data = guidata(gui_fig);

opts.Default = 'Yes';
opts.Interpreter = 'tex';
user_confirm = questdlg('\fontsize{14} Save?','Confirm exit',opts);
switch user_confirm
    case 'Yes'
        % Save and close
        atlas2histology_tform = ...
            gui_data.histology_ccf_manual_alignment;
        histology_control_points = gui_data.histology_control_points;
        atlas_control_points     = gui_data.atlas_control_points;
        save_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');

        save(save_fn,'atlas2histology_tform', 'atlas_control_points', 'histology_control_points');
        disp(['Saved ' save_fn]);
        delete(gui_fig);

    case 'No'
        % Close without saving
        delete(gui_fig);

    case 'Cancel'
        % Do nothing

end   

% Update toolbar GUI
% cp_lightsheet.update_toolbar_gui(gui_data.histology_toolbar_gui);

end


function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

idim    = gui_data.chooselist(gui_data.curr_slice, 2);
sluse   = gui_data.atlas_slice;
alldims = 1:3;
toplot  = find(alldims~=idim);
% 
% atlasize = size(gui_data.tv);
% if atlasize(idim) > sluse
%     sluse = atlasize(idim);
% end
% 
curr_atlas = volumeIdtoImage(gui_data.tv, [sluse idim]);

set(gui_data.atlas_im_h,'CData', curr_atlas);
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Reset histology-aligned atlas boundaries if not
% histology_aligned_atlas_boundaries_init = zeros(size(curr_image));
% set(gui_data.histology_aligned_atlas_boundaries, ...
%     'CData',histology_aligned_atlas_boundaries_init, ...
%     'AlphaData',histology_aligned_atlas_boundaries_init);

% Upload gui_data
guidata(gui_fig, gui_data);

end
















