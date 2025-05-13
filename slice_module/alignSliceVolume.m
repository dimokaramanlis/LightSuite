function [slicevol, sliceinfo] = alignSliceVolume(slicevol, sliceinfo)
%ALIGNSLICEVOLUME Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
if ~isnumeric(slicevol)
    % it is a path and we have to load it as a path
    assert(isstring(slicevol) | ischar(slicevol))
    slicevol = loadLargeSliceVolume(slicevol);
end
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% we first reorder the data volume
orderfile = fullfile(sliceinfo.procpath, 'volume_for_ordering_processing_decisions.txt');
if exist(orderfile, "file")
    tabledecisions = readtable(orderfile);
    sliceorder     = tabledecisions.NewOrderOriginalIndex;
    flipsdo        = logical(tabledecisions.FlipState);
else
    sliceorder = 1:sliceinfo.Nslices;
    flipsdo    = false(size(sliceorder));
end
slicevol(:, :, :, flipsdo) = flip(slicevol(:, :, :, flipsdo), 2);
slicevol                   = slicevol(:, :, :, sliceorder);

Nchans   = size(slicevol, 3);
%--------------------------------------------------------------------------
fprintf('Standardizing and filtering volume... '); tic;
% we standarding the volume values
ireg        = 1;
regdown     = 4;
scalefilter = 100/sliceinfo.px_process;
volregister = squeeze(slicevol(:, :, ireg, :));
sizedown    = [ceil(sliceinfo.size_proc/regdown) sliceinfo.Nslices];
volregister = imresize3(volregister, sizedown);

if sliceinfo.use_gpu
    volregister = single(gpuArray(volregister));
else
    volregister = single(volregister);
end
bval        = median(single(sliceinfo.backvalues(ireg,:)));
volregister = (volregister - bval)/bval;
volregister(volregister<0) = 0;

% we perform spatial bandpass filtering
imhigh       = spatial_bandpass(volregister, scalefilter/regdown, 3, 3, sliceinfo.use_gpu);
thresuse     = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp        = find(imhigh>thresuse);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Calculating point clouds... '); tic;
% we find all points above threshold
[row,col,slice] = ind2sub(size(imhigh), idxcp);
[~, ~, ic]      = unique(slice);
Nregister       = 1500;
rng(1);
allclouds       = cell(sliceinfo.Nslices, 1); % second dimension is flipped
sigmause        = sliceinfo.slicethickness/sliceinfo.px_process/4;

for ii = 1:sliceinfo.Nslices
    idxcurr     = ic == ii; 
    currx       = col(idxcurr)*regdown;
    curry       = row(idxcurr)*regdown;
    if sliceinfo.use_gpu
        currx       = gather(currx);
        curry       = gather(curry);
    end

    randz       = randn(nnz(idxcurr), 1)*sigmause;

    Npts      = nnz(idxcurr);
    randomfac = Nregister/Npts;

    pccurr        = pointCloud([currx curry randz]);
    allclouds{ii} = pcdownsample(pccurr,     'random', randomfac, 'PreserveStructure', true);
end
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% we here calculate alignments between slices
itemplate = ceil(sliceinfo.Nslices/2);
fprintf('==================Alignment forward pass==================\n')
indsforward  = itemplate:sliceinfo.Nslices;
transformsforward = alignConsecutiveSlices(allclouds, indsforward);
fprintf('==================Alignment backward pass==================\n')
indsbackwards = itemplate:-1:1;
transformsbackward = alignConsecutiveSlices(allclouds, indsbackwards);
%--------------------------------------------------------------------------
fprintf('Applying rigid transformations to generate aligned volume... '); tic;
% we then apply the transforms to the volume
slicevol(:, :, :, indsforward)  = applyConsecutiveTransforms(slicevol, ...
    indsforward, transformsforward, median(sliceinfo.backvalues,2));
slicevol(:, :, :, indsbackwards) = applyConsecutiveTransforms(slicevol, ...
    indsbackwards, transformsbackward, median(sliceinfo.backvalues,2));
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Saving aligned volume... '); tic;
saveLargeSliceVolume(slicevol, sliceinfo.channames, sliceinfo.slicevolfin);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Generating downsampled volume and saving... '); tic;

dpsavelowres = fullfile(sliceinfo.procpath, 'volume_for_inspection.tiff');
scalesize    = [ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register) sliceinfo.Nslices];
voldown      = zeros([scalesize(1:2) 3 scalesize(3)], 'uint8');
chansmap     = [3 2 1];

for ichan = 1:min(Nchans, 3)
    volproc     = imresize3(squeeze(slicevol(:, :, ichan, :)), scalesize);
    backproc    = single(median(sliceinfo.backvalues(ichan, :)));
    volproc     = (single(volproc) - backproc)./backproc;
    maxval      = quantile(volproc, 0.999, 'all');
    voldown(:, :, chansmap(ichan), :) =  uint8(255*volproc/maxval);
end

options.compress = 'lzw';
options.message  = false;
options.color    = true;
options.big      = false;

if exist(dpsavelowres, 'file')
    delete(dpsavelowres);
end
saveastiff(voldown, dpsavelowres, options);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
end

