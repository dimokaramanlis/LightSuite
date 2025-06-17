function [slicevol, regopts] = alignSliceVolume(slicevol, sliceinfo)
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
howtoperm   = [3 1 2];
ireg        = 1;
volregister = squeeze(slicevol(:, :, ireg, :));
scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
sizedown    = [finsize sliceinfo.Nslices ];
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
imhigh            = spatial_bandpass(volregister, scalefilter, 3, 3, sliceinfo.use_gpu);
imhigh            = permute(imhigh, howtoperm);
thresuse          = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
% thresuse          = min(6, thresuse);
idxcp             = find(imhigh>thresuse);
[slice, row, col] = ind2sub(size(imhigh), idxcp);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Loading and processing Allen Atlas template... '); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tv               = tv(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
av               = av(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
tvreg            = imresize3(tv, sliceinfo.px_atlas/sliceinfo.px_register);
avreg            = imresize3(av, sliceinfo.px_atlas/sliceinfo.px_register);
tv_points        = extractHighSFVolumePoints(tvreg, sliceinfo.px_register);
tv_cloud         = pointCloud(tv_points);
% Npts             = size(tv_points, 1);
% tv_cloud_use     = pcdownsample(pointCloud(tv_points),'random', 20000/Npts, 'PreserveStructure',true);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% we first align 3d volume to Allen

atlasframe = size(tv, [2 3]);
yvals      = slice * sliceinfo.slicethickness/sliceinfo.px_register;
xvals      = row;
zvals      = col;
Xinput     = [gather(xvals), gather(yvals), gather(zvals)];
pcsample   = pointCloud(Xinput);

errall = nan(3, 1);
tformslices(sliceinfo.Nslices, 1) = rigidtform2d;

for ii = 1:3
    if ii > 1
        tformslices     = refineSampleFromAtlas(tv_cloud, pcsample, tformrigid, 'rigid');
    end
    [tformrigid, errall(ii)] = alignAtlasToSample(tv_cloud, pcsample, tformslices);
end

slicevol = getRigidlyAlignedVolume(sliceinfo, slicevol, tformslices, atlasframe);

regopts = struct();
regopts.howtoperm                     = howtoperm;
regopts.tformrigid_allen_to_samp_20um = tformrigid;
regopts.howtoperm    = [3 1 2];
regopts.procpath     = sliceinfo.procpath;
regopts.registres    = sliceinfo.px_register;
regopts.allenres     = 10; % um
regopts.atlasaplims  = sliceinfo.atlasaplims;
regopts.pxsizes      = [sliceinfo.slicethickness/sliceinfo.px_register 1 1];
regopts.extentfactor = 6; % # slices to extend beyond rigid registration
save(fullfile(sliceinfo.procpath, 'regopts.mat'), '-struct', 'regopts')
%--------------------------------------------------------------------------
scalesize = [ceil(size(slicevol,[1 2])*sliceinfo.px_process/sliceinfo.px_register) sliceinfo.Nslices];
volsave   = single(gather(imresize3(squeeze(slicevol(:,:,1,:)), scalesize)));
regvolfac = (2^16-1)/max(volsave, [],"all");
voldown   = uint16(regvolfac*volsave);

samplepath = fullfile(sliceinfo.procpath, sprintf('sample_register_%dum.tif', regopts.registres));
options.compress = 'lzw';
options.message  = false;
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(voldown, samplepath, options);
%--------------------------------------------------------------------------
volsamp = single(permute(volregister, howtoperm))/single(quantile(volregister, 0.999,'all'));
rout    = imref3d([sliceinfo.Nslices size(avreg, [2 3])], 1, regopts.pxsizes(1), 1);
avtest  = imwarp(avreg, imref3d(size(avreg)), tformrigid, 'nearest', 'OutputView', rout);
for idim = 1:3
    cf = plotAnnotationComparison(uint8(255*volsamp), avtest, idim, regopts.pxsizes);
    savepngFast(cf, sliceinfo.procpath, ...
        sprintf('%s_dim%d_initial_registration', sliceinfo.mousename, idim), 300, 2);
    close(cf);
end
%--------------------------------------------------------------------------
fprintf('Saving aligned volume... '); tic;
saveLargeSliceVolume(slicevol, sliceinfo.channames, sliceinfo.slicevolfin);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Generating downsampled volume and saving... '); tic;

dpsavelowres = fullfile(sliceinfo.procpath, 'volume_for_inspection.tiff');
scalesize    = [ceil(size(slicevol,[1 2])*sliceinfo.px_process/sliceinfo.px_register) sliceinfo.Nslices];
voldown      = zeros([scalesize(1:2) 3 scalesize(3)], 'uint8');

for ichan = 1:min(Nchans, 3)
    volproc     = imresize3(squeeze(slicevol(:, :, ichan, :)), scalesize);
    backproc    = single(median(sliceinfo.backvalues(ichan, :)));
    volproc     = (single(volproc) - backproc)./backproc;
    maxval      = quantile(volproc, 0.999, 'all');
    voldown(:, :, ichan, :) =  uint8(255*volproc/maxval);
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

