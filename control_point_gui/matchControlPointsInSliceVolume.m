function matchControlPointsInSliceVolume(opts)
% Part of AP_histology toolbox
%
% Manually align histology slices and matched CCF slices

% STILL TO FIX: new predictions are not based on interpolation only

% Initialize guidata
gui_data = struct;
opts.downfac_reg = opts.allenres/opts.registres;
gui_data.save_path = opts.procpath;

%--------------------------------------------------------------------------
% Load atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
factv       = 255/single(max(gui_data.tv,[],"all"));
gui_data.tv = uint8(single(gui_data.tv(opts.atlasaplims(1):opts.atlasaplims(2), :, :))*factv);
gui_data.tv = imresize3(gui_data.tv,opts.downfac_reg);
gui_data.av = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
gui_data.av = imresize3(gui_data.av(opts.atlasaplims(1):opts.atlasaplims(2), :, :),...
    opts.downfac_reg, "Method","nearest");
gui_data.Rmoving  = imref3d(size(gui_data.av));
disp('Done.')
%--------------------------------------------------------------------------

% volume_dir       = dir(fullfile(opts.procpath,'sample_register_*um.tif'));
volume_dir       = dir(fullfile(opts.procpath,'*inspection.tif*'));

volpath          = fullfile(volume_dir.folder, volume_dir.name);
volload          = readDownStack(volpath);
volload          = permute(volload, [1 2 4 3]);
volload          = permute(volload, [opts.howtoperm 4]);
volload          = single(volload);
minvals          = single(quantile(volload, 0.01, [2 3]));
maxvals          = single(quantile(volload, 0.999, [2 3]));
volload          = 255 * (volload - minvals)./(maxvals - minvals);
gui_data.volume  = uint8(volload);

factv               = 255/single(max(gui_data.volume,[],"all"));
gui_data.volume     = uint8(single(gui_data.volume)*factv);
gui_data.Nslices    = size(gui_data.volume, 1);

chooselist          = generate_cp_list_slices(gui_data.volume);
gui_data.chooselist = chooselist;

% we start by placing the volume in the middle of the atlas
gui_data.slicepos   = opts.pxsizes(1)*(1:gui_data.Nslices)';
gui_data.slicepos   = gui_data.slicepos - median(gui_data.slicepos) + size(gui_data.tv, 1)/2;
gui_data.Ratlas     = imref3d(size(gui_data.tv));
gui_data.globalerr  = nan;
gui_data.NpointsTOT = 0;
gui_data.colsuse    = 1:size(volload, 4);
gui_data.avgspacing = nan;
%--------------------------------------------------------------------------
[rows, columns] = size(gui_data.tv, [2 3]);
linespacing     = round(min([rows, columns]/12));

xspaces = 1:linespacing:columns;
xlinesx = [repmat(xspaces, [2 1]); nan(1, numel(xspaces))];
xlinesy = [repmat([1; rows], [1 numel(xspaces)]) ; nan(1, numel(xspaces))];

yspaces    = 1:linespacing:rows;
ylinesy   = [repmat(yspaces, [2 1]); nan(1, numel(yspaces))];
ylinesx   = [repmat([1; columns], [1 numel(yspaces)]) ; nan(1, numel(yspaces))];
xlinesall = [xlinesx(:); ylinesx(:)];
ylinesall = [xlinesy(:); ylinesy(:)];
%--------------------------------------------------------------------------
% Load automated alignment
auto_ccf_alignment_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
if exist(auto_ccf_alignment_fn,'file')
    oldtform = load(auto_ccf_alignment_fn);
    % gui_data.histology_ccf_auto_alignment = oldtform.atlas2histology_tform;
    gui_data.histology_control_points     = oldtform.histology_control_points;
    gui_data.atlas_control_points         = oldtform.atlas_control_points;
else
    % Initialize alignment control points and tform matricies
    gui_data.histology_control_points = repmat({zeros(0,3)},size(gui_data.chooselist, 1),1);
    gui_data.atlas_control_points     = repmat({zeros(0,3)},size(gui_data.chooselist, 1),1);
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

gui_data.curr_slice = 1;
curr_image = squeeze(gui_data.volume(chooselist(gui_data.curr_slice, 1),:,:,:));
curr_image = blankImage_slice(curr_image, chooselist(gui_data.curr_slice, 2:end), true);
% curr_image = adapthisteq(curr_image(:, :, 1));

% Set up axis for histology image
gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
set(gui_data.histology_ax,'Position',[0,0,0.5,0.9]);
hold on; colormap(gray); axis image off;
gui_data.histology_im_h = image(curr_image,...
    'Parent',gui_data.histology_ax,'ButtonDownFcn',@mouseclick_histology);
% setup grid
gui_data.histology_grid = line(xlinesall(:), ylinesall(:), 'Color', 'w', 'LineWidth', 0.5);
set(gui_data.histology_grid,'Visible','off')

histology_aligned_atlas_boundaries_init = ...
    zeros(size(curr_image));
gui_data.histology_aligned_atlas_boundaries = ...
    plot(histology_aligned_atlas_boundaries_init(:,1), histology_aligned_atlas_boundaries_init(:,2),...
    'r.','MarkerSize',3, 'Parent',gui_data.histology_ax, 'PickableParts','none');

gui_data.atlas_slice = round(gui_data.slicepos(gui_data.curr_slice));
curr_atlas = squeeze(gui_data.tv(gui_data.atlas_slice, :, :));

% Set up axis for atlas slice
gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse'); 
set(gui_data.atlas_ax,'Position',[0.5,0,0.5,0.9]);
hold on; axis image off; colormap(gray); clim([0,255]);
gui_data.atlas_im_h = imagesc(curr_atlas, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);
% setup grid
gui_data.atlas_grid = line(xlinesall(:), ylinesall(:), 'Color', 'w', 'LineWidth', 0.5);
set(gui_data.atlas_grid,'Visible','off')

gui_data.histology_control_points_plot = plot(gui_data.histology_ax,nan,nan,...
    'o','MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'g');
gui_data.atlas_control_points_plot = plot(gui_data.atlas_ax,nan,nan,...
    'o','MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r');

% If there was previously auto-alignment, intitialize with that
if isfield(gui_data,'histology_ccf_auto_alignment')
    gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;
end

% Upload gui data
guidata(gui_fig,gui_data);

% Initialize alignment - FIX!!!
align_ccf_to_histology(gui_fig);

set_histology_title(gui_fig);

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Controls: \rm' ...
    'Left/right: switch slice' ...
    'click: set reference points for manual alignment (3 minimum)', ...
    'space: toggle alignment overlay visibility', ...
    'enter: go to slice id',...
    '1-3 or 0: toggle channels (0 for all)',...
    'g: toggle alignment grid', ...
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
    case 'g'
        curr_visibility = ...
            get(gui_data.histology_grid,'Visible');
        set(gui_data.histology_grid,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
        set(gui_data.atlas_grid,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
    
    case {'return', 'enter'}
        prompt = {sprintf('Move to slice (1-%d):', gui_data.Nslices)};
        dlg_title = 'Move to Position';
        answer = inputdlg(prompt, dlg_title, [1 60], {num2str(1)});
    
        if isempty(answer), return; end 
        new_pos = round(str2double(answer{1}));
    
        if isnan(new_pos) || new_pos < 1 || new_pos > gui_data.Nslices
            warndlg('Invalid input. Please enter a valid position number.', 'Input Error');
            return;
        end 

        slicego = find(gui_data.chooselist(:,1) == new_pos, 1);
        gui_data.curr_slice =slicego;
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % c: clear current points
    case 'c'
        gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,3);
        gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,3);
        
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
    case '1'
        gui_data.colsuse = 1;
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
    case '2'
        gui_data.colsuse = 2;
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
    case '3'
        gui_data.colsuse = 3;
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
    case '0'
        gui_data.colsuse = [1 2 3];
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % s: save
    case 's' 
        atlas2histology_tform = gui_data.histology_ccf_manual_alignment;
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


cpt     = zeros(1,3);
cpt([2 3]) = flip(eventdata.IntersectionPoint(1:2));
cpt(1)     = gui_data.chooselist(gui_data.curr_slice, 1);
toplot     = [2 3];


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

dval = gui_data.atlas_slice;
cpt     = zeros(1,3);
cpt([2 3]) =  flip(eventdata.IntersectionPoint(1:2));
cpt(1) = dval;

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    cpt);

toplot = [2 3];
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


Nmin = 4;
cptsatlas     = cat(1, gui_data.atlas_control_points{:});
cptshistology = cat(1, gui_data.histology_control_points{:});
sliceinds     = cptshistology(:, 1);
useglobaldata = false;
if size(cptshistology,1) == size(cptsatlas,1) && ...
        (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
    
    valsout = accumarray(cptshistology(:,1), cptsatlas(:,1), [gui_data.Nslices 1], @mean);
    iuse    = valsout>0;
    if nnz(iuse) > 1
        slicepos = interp1(find(iuse),valsout(iuse), 1:gui_data.Nslices, "linear", 'extrap');
        slicepos(iuse) = valsout(iuse);
        gui_data.slicepos = slicepos;
        gui_data.slicepos = slicepos;
    end
    useglobaldata = true;
end


cptsatlas     = cptsatlas(:, [2 1 3]);
cptshistology = cptshistology(:, [2 1 3]);

% we have to re-estimate the ap sampling based on user selection (could be
% non-uniform)

if size(cptshistology,1) == size(cptsatlas,1) && ...
        (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
    
    cptshistology(:, 2) = cptsatlas(:, 2); % there is no inherent spacing in histology...

    % If same number of >= 3 control points, use control point alignment
    [tform, mse] = fitRigidTrans3D(cptsatlas, cptshistology);
    % [tform,iidx] = estgeotform3d(cptsatlas, cptshistology, 'rigid', 'MaxDistance',10, 'MaxNumTrials',2e3);
    gui_data.volwrap = imwarp(gui_data.av, gui_data.Ratlas, tform, 'nearest','OutputView',gui_data.Ratlas);

    gui_data.globalerr  = mse;
    gui_data.NpointsTOT = size(cptsatlas, 1);

elseif isfield(gui_data,'histology_ccf_auto_alignment')
    % If no points, use automated outline if available
    tform = gui_data.histology_ccf_auto_alignment;
    tform = rigidtform3d(tform);
    gui_data.volwrap = imwarp(gui_data.av, gui_data.Rmoving, tform, 'nearest','OutputView',gui_data.Rmoving);
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
    tform = rigidtform3d;
    % allatlaspts = allatlaspts(:, [2 1 3]);
    % allpts = allpts(:, [2 1 3]);
    % tform = fitAffineTrans3D(allatlaspts, allpts);

    gui_data.volwrap = imwarp(gui_data.av,gui_data.Ratlas, tform, 'nearest', 'OutputView', gui_data.Ratlas);
end

atlasid            = round(gui_data.slicepos( gui_data.chooselist(gui_data.curr_slice,1)));
curr_slice_warp    = volumeIdtoImage(gui_data.volwrap, [atlasid 1]);

sliceid   =  gui_data.chooselist(gui_data.curr_slice,1);
idsforref = sliceinds == sliceid; 
if nnz(idsforref) > 3 & useglobaldata
    raim = imref2d(size(curr_slice_warp));
    histpts = cptshistology(idsforref, :);
    atpts   = cptsatlas(idsforref, :);
    newpts  = tform.transformPointsForward(atpts);
    taff  = fitgeotform2d(newpts(:, [3 1]), histpts(:, [3 1]),'affine');
    curr_slice_warp = imwarp(curr_slice_warp, raim, taff, 'nearest', 'Output', raim);
end

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

%--------------------------------------------------------------------------
% Get guidata
gui_data = guidata(gui_fig);
toplot     = [2 3];
%--------------------------------------------------------------------------
% Set next histology slice
idvol = gui_data.chooselist(gui_data.curr_slice, 1);
curr_image = squeeze(gui_data.volume(idvol, :, :, gui_data.colsuse));
[curr_image, iy, ix]= blankImage_slice(curr_image,  gui_data.chooselist(gui_data.curr_slice, 2:end), true);
% curr_image = adapthisteq(curr_image);

cptscurr = gui_data.atlas_control_points{gui_data.curr_slice};

if ~isempty(cptscurr)
    gui_data.atlas_slice = median(cptscurr(:,1));
else
    tform = rigidtform3d(gui_data.histology_ccf_manual_alignment);
    
    induse     = gui_data.chooselist(gui_data.curr_slice, 1);
    itform     = tform.invert;

    ieq = all(gui_data.chooselist(:,2:3) == gui_data.chooselist(gui_data.curr_slice, 2:3),2);
    atlaspred = cat(1, gui_data.atlas_control_points{ieq});
    histpred  = cat(1, gui_data.histology_control_points{ieq});
    unvals    = unique(histpred(:, 1));
    if numel(unvals) > 2
        valsout = accumarray(histpred(:,1), atlaspred(:,1), [gui_data.Nslices 1], @median);
        iuse = valsout>0;
        pfit  = interp1(find(iuse), valsout(iuse), induse, 'linear', 'extrap');
        appos = round(pfit);
    else
        appos      = round(gui_data.slicepos(induse));
    end
   


    
    [xl, yl, zl] = itform.outputLimits([min(iy) max(iy)], [appos appos], [min(ix) max(ix)]);
    sluse = round(median(yl));    
    sluse = max(sluse, 1);
    gui_data.atlas_slice = sluse;
end


% currlim    = getImageLimits(curr_image, 0.001);
set(gui_data.histology_im_h,'CData', curr_image)
% gui_data.histology_ax.CLim = [0 currlim(2)];
% Plot control points for slice
set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(1)));

histology_aligned_atlas_boundaries_init = nan(1,2);
set(gui_data.histology_aligned_atlas_boundaries, ...
    'XData',histology_aligned_atlas_boundaries_init(:,1), 'YData',histology_aligned_atlas_boundaries_init(:,2));


%--------------------------------------------------------------------------
% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

% update atlas slice
update_atlas_slice(gui_fig)

set_histology_title(gui_fig)
%--------------------------------------------------------------------------
end


function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point
gui_data.atlas_slice = ...
    gui_data.atlas_slice + eventdata.VerticalScrollCount;

gui_data.atlas_slice = max(gui_data.atlas_slice, 1);
gui_data.atlas_slice = min(gui_data.atlas_slice, size(gui_data.tv, 1));

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

end


function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);
%--------------------------------------------------------------------------
sluse   = gui_data.atlas_slice;
toplot  = [2 3];
 
% just in case something is wrong
atlasize = size(gui_data.tv, 1);
if atlasize < sluse
    sluse = atlasize;
end
%--------------------------------------------------------------------------
curr_atlas = volumeIdtoImage(gui_data.tv, [sluse 1]);
curr_atlas = adapthisteq(curr_atlas);

set(gui_data.atlas_im_h,'CData', curr_atlas);
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));
%--------------------------------------------------------------------------
% get an informative title
tform = rigidtform3d(gui_data.histology_ccf_manual_alignment);
[~, ~,~,tstr] = reportRotationAngles(tform.R, false);
tstr = sprintf('%s, NptsTOT = %d, mse = %2.1f', tstr, gui_data.NpointsTOT, gui_data.globalerr);
curraptxt = sprintf('Current AP atlas: %d um', sluse * 20);
title(gui_data.atlas_ax, {tstr, curraptxt});
%--------------------------------------------------------------------------
% Upload gui_data
guidata(gui_fig, gui_data);
%--------------------------------------------------------------------------
end


function set_histology_title(gui_fig)
gui_data = guidata(gui_fig);

currsliceid = gui_data.chooselist(gui_data.curr_slice, 1);
irelevant   = gui_data.chooselist(:, 1) == currsliceid;

ptscurr = cat(1, gui_data.histology_control_points{irelevant});

tstr = sprintf('Current slice: %d/%d, Npts = %d', currsliceid, gui_data.Nslices, size(ptscurr, 1));
title(gui_data.histology_ax, tstr);

end