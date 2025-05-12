function sliceinfo = getSliceInfo(filepaths)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


sliceinfo = struct();
Nfiles    = numel(filepaths);
allsizes  = cell(Nfiles, 2);

for ifile = 1:Nfiles
    dataim = BioformatsImage(filepaths{ifile});

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
sliceinfo.filepaths  = filepaths;

dataim               = BioformatsImage(filepaths{ifile});
dataim.series        = allsizes{ifile,2}(end);
sliceinfo.channames  = dataim.channelNames;

end