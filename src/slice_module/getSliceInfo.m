function sliceinfo = getSliceInfo(sliceinfo)
%GETSLICEINFO Populate sliceinfo fields by inspecting input files.
%
%   For CZI files:  pixel size and channel names are read automatically.
%   For TIF files:  sliceinfo.pxsize and sliceinfo.channames must be set
%                   by the caller before calling this function. Each TIF
%                   file in filepaths is treated as one slice (no scenes).
%                   Multi-channel TIFs should be multi-page (one page per
%                   channel); single-channel TIFs are also supported.

%--------------------------------------------------------------------------
dp        = fileparts(sliceinfo.filepaths{1});
procpath  = fullfile(dp, 'lightsuite'); makeNewDir(procpath);
%--------------------------------------------------------------------------
Nfiles    = numel(sliceinfo.filepaths);
allsizes  = cell(Nfiles, 2);

[~, ~, firstext] = fileparts(sliceinfo.filepaths{1});
is_tif = any(strcmpi(firstext, {'.tif', '.tiff'}));

if is_tif
    % ------------------------------------------------------------------
    % TIF mode: each file is one slice, dimensions read via imfinfo
    % ------------------------------------------------------------------
    if ~isfield(sliceinfo, 'pxsize') || isempty(sliceinfo.pxsize)
        error('For TIF input, sliceinfo.pxsize must be set manually (pixel size in um, e.g. [0.65 0.65]).');
    end

    fprintf('Reading list of .tif files... '); tiftic = tic;
    for ifile = 1:Nfiles
        info = imfinfo(sliceinfo.filepaths{ifile});
        allsizes{ifile, 1} = [info(1).Width, info(1).Height];
        allsizes{ifile, 2} = [1]; % one scene per TIF file
    end

    sliceinfo.maxsize   = flip(max(cat(1, allsizes{:,1})));
    sliceinfo.Nslices   = Nfiles;
    sliceinfo.sliceinds = allsizes;

    % Derive channel names from page count if not already set
    if ~isfield(sliceinfo, 'channames') || isempty(sliceinfo.channames)
        lastinfo = imfinfo(sliceinfo.filepaths{end});
        Npages   = numel(lastinfo);
        sliceinfo.channames = arrayfun(@(i) sprintf('ch%d', i), 1:Npages, 'UniformOutput', false);
    end

    fprintf('Done! Found %d valid slices. Took %2.1f s\n', sliceinfo.Nslices, toc(tiftic));

else
    % ------------------------------------------------------------------
    % CZI mode: use BioformatsImage to read scene metadata
    % ------------------------------------------------------------------
    fprintf('Reading list of .czi files... '); czifiletic = tic;
    for ifile = 1:Nfiles
        dataim = BioformatsImage(sliceinfo.filepaths{ifile});

        %select relevant series
        Nseries = dataim.seriesCount;
        allwh   = nan(Nseries, 2);
        allpx   = nan(Nseries, 2);
        for is = 1:Nseries
            dataim.series = is;
            allwh(is, :) = [dataim.width dataim.height];
            allpx(is, :) = dataim.pxSize;
        end
        irel = [1; find(diff(allwh(:,1))>0)+1];
        medw = median(allwh(irel, 1));
        irem = abs((allwh(irel,1) - medw)/medw)>0.3;
        irel(irem) = [];

        allsizes{ifile, 1} = allwh(irel, :);
        allsizes{ifile, 2} = irel;
    end

    sliceinfo.maxsize    = flip(max(cat(1, allsizes{:,1})));
    sliceinfo.Nslices    = sum(cellfun(@numel, allsizes(:,2)));
    sliceinfo.sliceinds  = allsizes;
    sliceinfo.pxsize     = dataim.pxSize;

    dataim               = BioformatsImage(sliceinfo.filepaths{ifile});
    dataim.series        = allsizes{ifile,2}(end);
    sliceinfo.channames  = dataim.channelNames;

    fprintf('Done! Found %d valid slices. Took %2.1f s\n', sliceinfo.Nslices, toc(czifiletic));
end

%--------------------------------------------------------------------------
sliceinfo.procpath    = procpath;
sliceinfo.volorder    = fullfile(sliceinfo.procpath, 'volume_for_ordering.tiff');
sliceinfo.slicevol    = fullfile(sliceinfo.procpath, 'volume_centered');
sliceinfo.slicevolfin = fullfile(sliceinfo.procpath, 'volume_aligned');
%--------------------------------------------------------------------------
end
