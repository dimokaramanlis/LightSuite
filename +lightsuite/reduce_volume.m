function reduce_volume(~,~,histology_toolbar_gui)
% Part of AP_histology toolbox
%
% Pad, center, and rotate slice images

% Get images (from path in GUI)
histology_toolbar_guidata = guidata(histology_toolbar_gui);

volume_dir = dir(fullfile(histology_toolbar_guidata.save_path,'*.tif'));
volpath    = fullfile(volume_dir.folder, volume_dir.name);
voldown    = readDownStack(volpath);

% Draw line to indicate midline for rotation
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7; % width/length
gui_width_fraction = 0.6; % fraction of screen width to occupy
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... % left x
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... % bottom y
    gui_width_px,gui_width_px/gui_aspect_ratio]; % width, height

gui_fig = figure('Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position);

align_axis = nan(3,3,2);
dids = 1:3;
for idim = 1:3
    currim = squeeze(quantile(voldown,0.9,idim));
    image(currim);
    colormap(gray)
    axis equal off;
    title('Draw a bounding rectangle (double-click when done)')
    curr_rect = drawrectangle;
    wait(curr_rect)
    currpos = curr_rect.Position; 
    rstart  = max(floor(currpos(1:2)), 1);
    rend    = min(ceil(currpos(1:2) + currpos(3:4)), flip(size(currim)));
    align_axis(idim,  ~(dids==idim), :) = [flip(rstart)' flip(rend)'];
end
close(gui_fig);

allends = round(squeeze(median(align_axis,1, 'omitnan')));
volred  = voldown(allends(1,1):allends(1,2), allends(2,1):allends(2,2), allends(3,1):allends(3,2));

% Make re-ordered filenames (with '_reorder to avoid overwriting)
reorder_target = strrep(volpath,'.tif','_reduced.tif');
saveLightsheetVolume(volred, reorder_target)
movefile(reorder_target, volpath);
% TO DO
% SAVE COORDINATES!!!
disp('Saved reduced volume');





