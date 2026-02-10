function  generateReordedVolume(sliceinfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

orderfile = fullfile(sliceinfo.procpath, 'volume_for_ordering_processing_decisions.txt');
if exist(orderfile, "file")
    tabledecisions = readtable(orderfile);
    sliceorder     = tabledecisions.NewOrderOriginalIndex;
    flipsdo        = tabledecisions.FlipState == 1;
    toremove       = tabledecisions.FlipState == -1;
else
    sliceorder = 1:sliceinfo.Nslices;
    flipsdo    = false(size(sliceorder));
    toremove   = false(size(sliceorder));
end
volread                = readDownStack(sliceinfo.volorder);
volread(:, :, :, flipsdo) = flip(volread(:, :, :, flipsdo), 2);
volread                = volread(:, :, :, sliceorder);
volread(:, :, :, toremove(sliceorder)) = []; % remove bad slices

options.compress = 'lzw';
options.message  = false;
options.color    = true;
options.big      = false;
dpsaveorder      = fullfile(sliceinfo.procpath, 'volume_ordered.tiff');
if exist(dpsaveorder, 'file')
    delete(dpsaveorder);
end
saveastiff(volread, dpsaveorder, options);

end