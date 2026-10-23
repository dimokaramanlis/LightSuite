function dp = resolveSampleFolder(datafolderpath, mousename)
%RESOLVESAMPLEFOLDER Find the folder that holds one sample's raw images.
%
%   DP = resolveSampleFolder(DATAFOLDERPATH, MOUSENAME) supports two layouts:
%     1. DATAFOLDERPATH contains one subfolder per mouse, with MOUSENAME in
%        its name (e.g. D:\DATA\20250101_DK001\). That subfolder is returned.
%     2. DATAFOLDERPATH itself holds the images (.czi/.tif/.tiff). It is
%        returned as is.
%
%   Only folders are matched, so image files named after the mouse are
%   never mistaken for the sample folder. Errors if the match is ambiguous.

    datafolderpath = char(datafolderpath);
    if ~exist(datafolderpath, 'dir')
        error('resolveSampleFolder:notFound', 'Data folder not found: %s', datafolderpath);
    end

    cands = [];
    if ~isempty(mousename)
        cands = dir(fullfile(datafolderpath, sprintf('*%s*', mousename)));
        cands = cands([cands.isdir] & ~ismember({cands.name}, {'.', '..'}));
    end

    if numel(cands) == 1
        dp = fullfile(cands.folder, cands.name);
    elseif numel(cands) > 1
        % prefer an exact name match, otherwise ask the user to be specific
        iexact = strcmp({cands.name}, mousename);
        if nnz(iexact) == 1
            dp = fullfile(cands(iexact).folder, cands(iexact).name);
        else
            error('resolveSampleFolder:ambiguous', ...
                ['Several folders in %s match ''%s'':\n  %s\n' ...
                 'Make mousename more specific or point datafolderpath to the sample folder.'], ...
                datafolderpath, mousename, strjoin({cands.name}, '\n  '));
        end
    elseif hasImages(datafolderpath)
        dp = datafolderpath;
    else
        error('resolveSampleFolder:notFound', ...
            ['No subfolder matching ''%s'' and no .czi/.tif(f) images in %s.\n' ...
             'Set datafolderpath to the folder with the images, or to its parent folder.'], ...
            mousename, datafolderpath);
    end
    fprintf('Sample folder: %s\n', dp);
end

function tf = hasImages(folder)
    tf = ~isempty(dir(fullfile(folder, '*.czi'))) || ...
         ~isempty(dir(fullfile(folder, '*.tif'))) || ...
         ~isempty(dir(fullfile(folder, '*.tiff')));
end
