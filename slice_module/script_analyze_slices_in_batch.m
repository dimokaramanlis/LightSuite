
% folder which contains mouse subfolders
datafolderpath = 'D:\Histology\Experiments\TRAP\'; %'D:\example_charlie'; %'J:\'; % 

mice = dir(datafolderpath);
mice = mice([mice.isdir] & ~startsWith({mice.name},'.'));

for mouseNum = numel(mice):-1:1
    mousename      = mice(mouseNum).name;
    
    dp = fullfile(datafolderpath, sprintf('*%s*', mousename));
    dp = dir(dp);
    dp = fullfile(dp.folder, dp.name);
    sliceinfo = parseSettingsFile(fullfile(dp, 'local_settings.txt'));
    
    sliceinfo.mousename   = mousename;
    filelistcheck         = dir(fullfile(dp, '*.czi'));
    filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
    sliceinfo.filepaths   = filepaths;
    sliceinfo             = getSliceInfo(sliceinfo);
    
    %% (auto) we first generate the slice volume
    % slicevol = generateSliceVolume(sliceinfo, sliceinfo.regchan);
      %% (manual) reorder, flip and discard slices if needed
    % SliceOrderEditor(sliceinfo.volorder)
    generateReordedVolume(sliceinfo);
    % % (auto) we align slices and initialize registration
    % sliceinfo          = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
    % sliceinfo          = sliceinfo.sliceinfo;
    % sliceinfo          = copyStructBtoA(sliceinfo, settings);
    % alignedvol         = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

    %% Determine angle
    opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
    determineCuttingAngleGUI(opts)
end

% end

%     %% (auto) detect cells in slices
%     sliceinfo  = load(fullfile(sliceinfo.procpath, "sliceinfo.mat"));
%     sliceinfo  = sliceinfo.sliceinfo;
%     ichan      = find(contains(sliceinfo.channames, 'Cy3')); % We use the tdTomato channel!
%     sliceinfo.debug    = true; % toggle plotting (takes longer) for detections
%     sliceinfo.celldiam = 14; % expected cell diameter in um
%     sliceinfo.thresuse = single([0.75 0.4]); % thresholds in SBR units (first for detection and then for expansion)
%     celllocs = extractCellsFromSliceVolume(sliceinfo, ichan);
% 
%     fsavename        = fullfile(opts.savepath, 'cell_locations_slice.mat');
%     save(fsavename, 'celllocs');
% end