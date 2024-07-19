function cp_lightsheet
% Toolbar GUI for running histology pipeline

% Set up the gui
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.5; % width/length
gui_width_fraction = 0.4; % fraction of screen width to occupy
gui_border = 50; % border from gui to screen edge
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_height_px = gui_width_px/gui_aspect_ratio;
gui_position = [...
    gui_border, ... % left x
    screen_size_px(4)-(gui_height_px+gui_border+50), ... % bottom y
    gui_width_px,gui_height_px]; % width, height

histology_toolbar_gui = figure('Toolbar','none','Menubar','none','color','w', ...
    'Name','Control-point registration for lightsheet', ...
    'Units','pixels','Position',gui_position);

% Set up the text to display coordinates
gui_data.gui_text = annotation('textbox','String','','interpreter','tex', ...
    'Units','normalized','Position',[0,0,1,1],'VerticalAlignment','top', ...
    'FontSize',12,'FontName','Consolas','PickableParts','none');

% File menu
gui_data.menu.file = uimenu(histology_toolbar_gui,'Text','File selection');
uimenu(gui_data.menu.file,'Text','Set raw image path','MenuSelectedFcn',{@set_image_path,histology_toolbar_gui});
uimenu(gui_data.menu.file,'Text','Set processing save path','MenuSelectedFcn',{@set_save_path,histology_toolbar_gui});

% Preprocessing menu
gui_data.menu.preprocess = uimenu(histology_toolbar_gui,'Text','Volume preprocessing');
uimenu(gui_data.menu.preprocess,'Text','Create alignment volume','MenuSelectedFcn', ...
    {@cp_lightsheet.create_slice_images,histology_toolbar_gui});
uimenu(gui_data.menu.preprocess,'Text','Reduce volume','MenuSelectedFcn', ...
    {@cp_lightsheet.reduce_volume,histology_toolbar_gui});
uimenu(gui_data.menu.preprocess,'Text','Re-order dimensions','MenuSelectedFcn', ...
    {@cp_lightsheet.reorder_dimensions,histology_toolbar_gui});
% uimenu(gui_data.menu.preprocess,'Text','Orient volume','MenuSelectedFcn', ...
%     {@cp_lightsheet.orient_volume,histology_toolbar_gui});

% Atlas menu
gui_data.menu.atlas = uimenu(histology_toolbar_gui,'Text','Atlas alignment');
uimenu(gui_data.menu.atlas,'Text','Manually align volume (coronal)','MenuSelectedFcn', ...
    {@cp_lightsheet.register_volume,histology_toolbar_gui, 'coronal'});
uimenu(gui_data.menu.atlas,'Text','Manually align volume (axial)','MenuSelectedFcn', ...
    {@cp_lightsheet.register_volume,histology_toolbar_gui, 'axial'});
uimenu(gui_data.menu.atlas,'Text','Manually align volume (sagittal)','MenuSelectedFcn', ...
    {@cp_lightsheet.register_volume,histology_toolbar_gui, 'sagittal'});
uimenu(gui_data.menu.atlas,'Text','Add control points','MenuSelectedFcn', ...
    {@cp_lightsheet.add_control_points,histology_toolbar_gui});

% Annotation menu
gui_data.menu.annotation = uimenu(histology_toolbar_gui,'Text','Annotation');
uimenu(gui_data.menu.annotation,'Text','Neuropixels probes','MenuSelectedFcn', ...
    {@cp_lightsheet.annotate_neuropixels,histology_toolbar_gui});

% View menu
gui_data.menu.view = uimenu(histology_toolbar_gui,'Text','View');
uimenu(gui_data.menu.view,'Text','View aligned histology','MenuSelectedFcn', ...
    {@cp_lightsheet.view_aligned_histology,histology_toolbar_gui});

% Create GUI variables
gui_data.image_path = char;
gui_data.save_path = char;

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
cp_lightsheet.update_toolbar_gui(histology_toolbar_gui);

end

function set_image_path(h,eventdata,histology_toolbar_gui)

% Get guidata
gui_data = guidata(histology_toolbar_gui);

% Pick image path
% gui_data.image_path = uigetdir([],'Select path with raw images');
gui_data.image_path = 'D:\DATA_folder\Mice\DK031\Anatomy';

% Clear processed path (if there's one selected)
gui_data.save_path = [];

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
cp_lightsheet.update_toolbar_gui(histology_toolbar_gui);

end

function set_save_path(h,eventdata,histology_toolbar_gui)

% Get guidata
gui_data = guidata(histology_toolbar_gui);

% Pick image path
% gui_data.save_path = uigetdir([],'Select path to save processing');
gui_data.save_path = 'D:\DATA_folder\Mice\DK031\Anatomy';

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
cp_lightsheet.update_toolbar_gui(histology_toolbar_gui);

end




















