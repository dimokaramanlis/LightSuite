function runJetrawCompressColm(pathtarget)
    % Construct the base command
    %======================================================================\
    dplist = uigetdir('C:\', 'Select imaging folder');
    if ~exist(pathtarget, 'dir'), mkdir(pathtarget), end
    %======================================================================
    copyfilesFromFolder(dplist, pathtarget, '.ini');
    copyfilesFromFolder(dplist, pathtarget, '.txt');
    copyfilesFromFolder(dplist, pathtarget, '.xml');
    %======================================================================
    vw0folder = fullfile(pathtarget, 'VW0');
    if ~exist(vw0folder, 'dir'), mkdir(vw0folder), end
    %======================================================================   
    inidata    = inifile(fullfile(dplist, 'Experiment.ini'), 'readall');
    [row, ~] = find(contains(inidata, 'laser (nm)'));
    linedata   = inidata(row,:);
    chancands  = cellfun(@(x) str2double(x(end)), linedata(:,1));
    chanId     = chancands(contains(linedata(:, end), '561'));
    chanstr    = sprintf('*CHN%02d_*.tif', chanId);
    %======================================================================
    pathList           = dir(fullfile(dplist, '**',chanstr));
    [allfolders, ~, ~] = unique({pathList(:).folder});
    %======================================================================
    fprintf('Compressing tiffs...\n')
    Nfolders  = numel(allfolders);
    msg = []; tic;
    for ifolder = 1:Nfolders
        
        [~, finpath] = fileparts(allfolders{ifolder});
        folderpath   = fullfile(vw0folder, finpath);
        if ~exist(folderpath, 'dir'), mkdir(folderpath), end
        
        commandStr = sprintf('%s %s %s', 'jetraw compress -d',...
                folderpath, allfolders{ifolder});
        % Trim any trailing space
        commandStr = strtrim(commandStr);
        % Run the command
        [status, cmdout] = system(commandStr);

%         
%         currids     = iun==ifolder;
%         Nfiles      = nnz(currids);
%         batchsize   = min(1, Nfiles);
%         Nbatches    = ceil(Nfiles/batchsize);
%         folderfiles = filenames(currids);
%         
%         for ibatch = 1:Nbatches
%             istart = (ibatch - 1)*batchsize + 1;
%             iend   = min(ibatch * batchsize, Nfiles);
%             
%             commandStr = sprintf('%s %s %s', 'jetraw decompress -d',...
%                 folderpath, folderfiles{istart:iend});
%             % Trim any trailing space
%             commandStr = strtrim(commandStr);
%             % Run the command
%             [status, cmdout] = system(commandStr);
%         end
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Tile %d/%d. Time elapsed %2.2f s...\n', ifolder, Nfolders,toc);
        fprintf(msg);
    end
    %======================================================================


   

end


function copyfilesFromFolder(folderin, folderout, ext)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

folderstuff = dir(fullfile(folderin, ['*' ext]));
for ii = 1:numel(folderstuff)
    sourcefile = fullfile(folderstuff.folder, folderstuff.name);
    targetfile = fullfile(folderout, folderstuff.name);
    copyfile(sourcefile, targetfile)
end


end

