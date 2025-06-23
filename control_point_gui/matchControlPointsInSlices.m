function matchControlPointsInSlices(opts)
% Part of AP_histology toolbox
%
% Manually align histology slices and matched CCF slices

% Initialize guidata
gui_data = struct;
opts.downfac_reg = opts.allenres/opts.registres;

% Load atlas
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
factv            = 255/single(max(gui_data.tv,[],"all"));
gui_data.tv      = uint8(single(gui_data.tv(opts.atlasaplims(1):opts.atlasaplims(2), :, :))*factv);
gui_data.tv      = imresize3(gui_data.tv,opts.downfac_reg);
gui_data.av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
gui_data.av      = imresize3(gui_data.av(opts.atlasaplims(1):opts.atlasaplims(2), :, :),...
    opts.downfac_reg, "Method","nearest");
disp('Done.')

gui_data.save_path = opts.procpath;

volume_dir       = dir(fullfile(opts.procpath,'sample_register_*um.tif'));
volpath          = fullfile(volume_dir.folder, volume_dir.name);
volload          = readDownStack(volpath, 1);
volload          = permute(volload, opts.howtoperm);
gui_data.volume  = volload;

factv            = 255/single(max(gui_data.volume,[],"all"));
gui_data.volume  = uint8(single(gui_data.volume)*factv);
gui_data.Nslices = size(gui_data.volume, 1);

nfac    = opts.extentfactor;
Ratlas  = imref3d(size(gui_data.tv));
Rvolume = imref3d(size(gui_data.volume), 1, opts.pxsizes(1), 1);
yworld  = [Rvolume.YWorldLimits(1)-opts.pxsizes(1)*nfac, Rvolume.YWorldLimits(2)+nfac*opts.pxsizes(1)];
ypix    = ceil(range(yworld));
Rout    = imref3d([ypix, size(gui_data.tv, [2 3])], Rvolume.XWorldLimits, yworld,Rvolume.ZWorldLimits);
gui_data.tv = imwarp(gui_data.tv, Ratlas, opts.tformrigid_allen_to_samp_20um, 'linear',  'OutputView', Rout);
gui_data.av = imwarp(gui_data.av, Ratlas, opts.tformrigid_allen_to_samp_20um, 'nearest', 'OutputView', Rout);
% 
% Ratlas         = imref3d(size(gui_data.av));
% Rvolume        = imref3d(size(volload), 1, opts.pxsizes(1), 1);
% [avnew,    ~]  = imwarp(gui_data.av, Ratlas, opts.tformrigid_allen_to_samp_20um, 'nearest');
% [tvnew, rnew]  = imwarp(gui_data.tv, Ratlas, opts.tformrigid_allen_to_samp_20um, 'linear');
% 
% xworld     = [Rvolume.XWorldLimits(1) + 0.5, Rvolume.XWorldLimits(2) - 0.5];
% zworld     = [Rvolume.ZWorldLimits(1) + 0.5, Rvolume.ZWorldLimits(2) - 0.5];
% zworld(1)  = max(zworld(1), rnew.ZWorldLimits(1));
% zworld(2)  = min(zworld(2), rnew.ZWorldLimits(2));
% yworld     = [Rvolume.YWorldLimits(1) + 0.5 - 10*opts.pxsizes(1),...
%     Rvolume.YWorldLimits(2) + 10*opts.pxsizes(1) - 0.5];
% yworld(1)  = max(yworld(1), rnew.YWorldLimits(1));
% yworld(2)  = min(yworld(2), rnew.YWorldLimits(2));
% 
% [yy,xx,zz]     = rnew.worldToSubscript(xworld, yworld, zworld);
% gui_data.av    = avnew(yy(1):yy(2), xx(1):xx(2), zz(1):zz(2));
% gui_data.tv    = tvnew(yy(1):yy(2), xx(1):xx(2), zz(1):zz(2));
% gui_data.tv    = smoothdata(gui_data.tv, 1, 'gaussian', opts.pxsizes(1)/2);
yatlasvals     = linspace(yworld(1), yworld(2), ypix + 1);
yatlasvals     = yatlasvals(1:end-1) + median(diff(yatlasvals))/2;
ysamplevals    = linspace(Rvolume.YWorldLimits(1), Rvolume.YWorldLimits(2), gui_data.Nslices+1);
ysamplevals    = ysamplevals(1:end-1) + median(diff(ysamplevals))/2;
[~, atlasinds] = min(pdist2(ysamplevals',yatlasvals'), [],2);

gui_data.atlasinds  = atlasinds;
gui_data.slicewidth = median(diff(yatlasvals))* opts.pxsizes(1);
gui_data.atlasvals  = yatlasvals;
gui_data.Rmoving    = imref2d(size(gui_data.av, [2 3]));
gui_data.Rfixed     = imref2d(size(volload, [2 3]));


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
    % gui_data.histology_ccf_auto_alignment = oldtform.atlas2histology_tform;
    gui_data.histology_control_points     = oldtform.histology_control_points;
    gui_data.atlas_control_points         = oldtform.atlas_control_points;
else
    % Initialize alignment control points and tform matricies
    gui_data.histology_control_points = repmat({zeros(0,3)}, gui_data.Nslices, 1);
    gui_data.atlas_control_points     = repmat({zeros(0,3)}, gui_data.Nslices,1);
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
curr_image = volumeIdtoImage(gui_data.volume, [gui_data.curr_slice 1]);
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

gui_data.atlas_slice = gui_data.atlasinds(1);
curr_atlas = squeeze(gui_data.tv(gui_data.atlas_slice, :, :));

% Set up axis for atlas slice
gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse'); 
set(gui_data.atlas_ax,'Position',[0.5,0,0.5,0.9]);
hold on; axis image off; colormap(gray); clim([0,250]);
gui_data.atlas_im_h = imagesc(curr_atlas, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);

title(gui_data.atlas_ax, sprintf('Atlas slice = %2.2f h-slice widths', ...
        gui_data.atlas_slice/gui_data.slicewidth));

gui_data.histology_control_points_plot = plot(gui_data.histology_ax,nan,nan,'.g','MarkerSize',20);
gui_data.atlas_control_points_plot = plot(gui_data.atlas_ax,nan,nan,'.r','MarkerSize',20);

% If there was previously auto-alignment, intitialize with that
if isfield(gui_data,'histology_ccf_auto_alignment')
    gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;
end

% Upload gui data
guidata(gui_fig,gui_data);

% Initialize alignment - FIX!!!
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
            min(gui_data.curr_slice + 1, gui_data.Nslices);
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
        histology_control_points = gui_data.histology_control_points;
        atlas_control_points     = gui_data.atlas_control_points;
        save_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
        save(save_fn,'atlas_control_points', 'histology_control_points');
        disp(['Saved ' save_fn]);
        
end

end


function mouseclick_histology(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);
toplot   = [2 3];
cpt(1)   = gui_data.curr_slice;
cpt(2:3) = flip(eventdata.IntersectionPoint(1:2));


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
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology(gui_fig)
end

end


function mouseclick_atlas(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

cpt      = zeros(1,3);
cpt(1)   = gui_data.atlas_slice;
cpt(2:3) =  flip(eventdata.IntersectionPoint(1:2));
toplot   = [2 3];

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    cpt);

set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology(gui_fig)
end

end


function align_ccf_to_histology(gui_fig)

% Get guidata 
gui_data = guidata(gui_fig);


Nmin = 3;
cptsatlas     = gui_data.atlas_control_points{gui_data.curr_slice};
idatlas       = gui_data.atlas_slice;
cptsatlas     = cptsatlas(:, [3 2]); 
cptshistology = gui_data.histology_control_points{gui_data.curr_slice};
cptshistology = cptshistology(:, [3 2]);
currim        = squeeze(gui_data.av(idatlas, :, :));

tstrcurr      = sprintf('Slice %d/%d', gui_data.curr_slice, gui_data.Nslices);

if size(cptshistology,1) == size(cptsatlas,1) && ...
        (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
    
    tform = fitgeotform2d(cptsatlas, cptshistology, "affine");
    mse = mean(sqrt(sum((cptshistology-tform.transformPointsForward(cptsatlas)).^2, 2)));

    % tform = estgeotform3d(cptsatlas, cptshistology, 'similarity');
    currim = imwarp(currim, gui_data.Rmoving, tform, 'OutputView',gui_data.Rfixed);

    tstrcurr = sprintf('%s, Npoints = %d, mse = %2.2f ', tstrcurr, size(cptshistology,1), mse);

elseif size(gui_data.histology_control_points{gui_data.curr_slice},1) >= 1 ||  ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) >= 1
    % If less than 3 or nonmatching points, use auto but don't draw
    tstrcurr = sprintf('%s, New alignment ', tstrcurr);

    % Upload gui data
    guidata(gui_fig, gui_data);
    return

elseif isfield(gui_data,'histology_ccf_auto_alignment')
    % If no points, use automated outline if available
    tform = gui_data.histology_ccf_auto_alignment;
    tform = affinetform2d(tform);
    tstrcurr = sprintf('%s, Previous alignment ', tstrcurr);
    currim = imwarp(currim, gui_data.Rmoving, tform, 'OutputView',gui_data.Rfixed);
else
    tform = affinetform2d;
end

title(gui_data.histology_ax, tstrcurr)

av_warp_boundaries = round(conv2(currim,ones(3)./9,'same')) ~= currim;

[row,col] = ind2sub(size(currim), find(av_warp_boundaries));


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

cpointsatlas = gui_data.atlas_control_points{gui_data.curr_slice};
if ~isempty(cpointsatlas)
    gui_data.atlas_slice = median(cpointsatlas(:,1));
else
    gui_data.atlas_slice = gui_data.atlasinds(gui_data.curr_slice);
end


toplot     = [2 3];

curr_image = volumeIdtoImage(gui_data.volume, [gui_data.curr_slice 1]);
curr_image = adapthisteq(curr_image);
currlim    = getImageLimits(curr_image, 0.001);

set(gui_data.histology_im_h,'CData', curr_image)
gui_data.histology_ax.CLim = [0 currlim(2)];
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

% Update toolbar GUI
% cp_lightsheet.update_toolbar_gui(gui_data.histology_toolbar_gui);

end


function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

toplot     = [2 3];
sluse      = gui_data.atlas_slice;
curr_atlas = squeeze(gui_data.tv(sluse, :, :));
curr_atlas = adapthisteq(curr_atlas);

set(gui_data.atlas_im_h,'CData', curr_atlas);
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Reset histology-aligned atlas boundaries if not
% histology_aligned_atlas_boundaries_init = zeros(size(curr_image));
% set(gui_data.histology_aligned_atlas_boundaries, ...
%     'CData',histology_aligned_atlas_boundaries_init, ...
%     'AlphaData',histology_aligned_atlas_boundaries_init);


title(gui_data.atlas_ax, sprintf('Atlas slice = %2.2f h-slice widths', ...
        sluse/gui_data.slicewidth));

% Upload gui_data
guidata(gui_fig, gui_data);

end
















