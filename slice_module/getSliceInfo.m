function sliceinfo = getSliceInfo(sliceinfo)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
dp        = fileparts(sliceinfo.filepaths{1});
procpath  = fullfile(dp, 'lightsuite'); makeNewDir(procpath);
%--------------------------------------------------------------------------
Nfiles    = numel(sliceinfo.filepaths);
allsizes  = cell(Nfiles, 2);
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
sliceinfo.procpath   = procpath;

sliceinfo.volorder    = fullfile(sliceinfo.procpath, 'volume_for_ordering.tiff');
sliceinfo.slicevol    = fullfile(sliceinfo.procpath, 'volume_centered');
sliceinfo.slicevolfin = fullfile(sliceinfo.procpath, 'volume_aligned');
fprintf('Done! Found %d valid slices. Took %2.1f s\n', sliceinfo.Nslices, toc(czifiletic)); 
%--------------------------------------------------------------------------
end