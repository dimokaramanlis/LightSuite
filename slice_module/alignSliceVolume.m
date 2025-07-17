function [slicevol, regopts] = alignSliceVolume(slicevol, sliceinfo)
%ALIGNSLICEVOLUME Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
orderfile = fullfile(sliceinfo.procpath, 'volume_for_ordering_processing_decisions.txt');
if ~isnumeric(slicevol)
    % it is a path and we have to load it as a path
    assert(isstring(slicevol) | ischar(slicevol))
    slicevol = loadLargeSliceVolume(slicevol, 1:numel(sliceinfo.channames), orderfile);
end
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Applying ordering/flipping... '); tic;
% we first reorder the data volume
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
slicevol(:, :, :, flipsdo) = flip(slicevol(:, :, :, flipsdo), 2);
slicevol                   = slicevol(:, :, :, sliceorder);
slicevol(:, :, :, toremove(sliceorder)) = []; % remove bad slices

[Nchans, Nslices]   = size(slicevol, [3,4]);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Standardizing and filtering volume... '); tic;
% we standarding the volume values
howtoperm   = [3 1 2];
ireg        = 1;
volregister = squeeze(slicevol(:, :, ireg, :));
scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
sizedown    = [finsize Nslices ];
medpx       = 2*floor((21./sliceinfo.px_process)/2) + 1;
sliceinfo.medianfiltreg  = getOr( sliceinfo, 'medianfiltreg', false);
sliceinfo.use_gpu  = getOr( sliceinfo, 'use_gpu', true);

if sliceinfo.medianfiltreg 
    for islice = 1:Nslices
        volregister(:, :, islice) = medfilt2(volregister(:, :, islice), medpx*[1 1]);
    end
end

volregister = imresize3(volregister, sizedown);

if sliceinfo.use_gpu
    volregister = single(gpuArray(volregister));
else
    volregister = single(volregister);
end
fprintf('Done! Took %2.2f s\n', toc); 
%==========================================================================
fprintf('Extracting point clouds from data... '); tic;
% volflat           = volregister;
% for ii = 1:Nslices
%     volflat(:,:,ii) = imflatfield(volregister(:,:,ii), 100);
% end
% we perform spatial bandpass filtering
imhigh            = spatial_bandpass(volregister, scalefilter, 6, 3, sliceinfo.use_gpu);
imhigh            = permute(imhigh, howtoperm);
thresuse          = quantile(imhigh,0.99, [2 3])/2;
ptsall            = cell(Nslices, 1);
Nrows = size(imhigh, 2);
Ncols = size(imhigh, 3);
pthres = 0.15;
for ii = 1:Nslices
    % find relevant points
    idxcp      = find(squeeze(imhigh(ii, :, :))>thresuse(ii));
    [row, col] = ind2sub(size(imhigh, [2 3]), idxcp);
    % filter points that are in straight lines (remove tiling artifacts)
    straighty = accumarray(row, 1, [Nrows 1], @sum);
    iremrow   = find(straighty/Ncols > pthres);
    straightx = accumarray(col, 1, [Ncols 1], @sum);
    iremcol   = find(straightx/Nrows > pthres);
    ikeep     = ~ismembc(col, iremcol) & ~ismembc(row, iremrow);

    ptsall{ii} = [ones(nnz(ikeep), 1)*ii row(ikeep), col(ikeep)];

    % imagesc(squeeze(imhigh(ii, :, :)), [thresuse(ii)*0.95 thresuse(ii)]); hold on;
    % plot(col(ikeep), row(ikeep), 'r.')
    % pause;
end

Xcloud = cat(1, ptsall{:});
Xcloud(:, 1) = Xcloud(:, 1)* sliceinfo.slicethickness/sliceinfo.px_register;
pcsample   = pointCloud(Xcloud(:,[2 1 3]));

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
tv_points        = extractHighSFVolumePoints(tvreg, sliceinfo.px_register, sliceinfo.use_gpu);
tv_cloud         = pointCloud(tv_points);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% optimization steps
tformslices(Nslices, 1) = rigidtform2d;
Tvec       = median(pcsample.Location)- range(tv_cloud.Location)/2;
tformrigid = rigidtform3d(eye(3),Tvec);

% tvclouddown = pcdownsample(tv_cloud, 'random', (2e4/tv_cloud.Count));
% pcclouddown = pcdownsample(pcsample, 'random', (2e4/pcsample.Count));
% 
% [tformrigid, movreg] = pcregisterndt(tvclouddown,pcclouddown,50, 'Verbose',true);

transformtypes = {'rigid', 'similarity', 'affine'};
transformsteps = [1 1 1];
errall = nan(numel(transformsteps), 1);
for istep = 1:numel(transformsteps)
    currtranstype               = transformtypes{transformsteps(istep)};
    fprintf('Optimization step %d/%d: %s\n', istep, numel(transformsteps), currtranstype)
    tformslices                 = refineSampleFromAtlas(tv_cloud, pcsample, tformrigid, currtranstype);
    [tformrigid, errall(istep)] = alignAtlasToSample(tv_cloud, pcsample, tformslices);
    [rrx, rry, rrz]             = reportRotationAngles(tformrigid.R);
    fprintf('%s\n', repmat('=', [1 75]));
end
%--------------------------------------------------------------------------
atlasframe = size(tv, [2 3]);
slicevol   = getRigidlyAlignedVolume(sliceinfo, slicevol, tformslices, atlasframe);

regopts = struct();
regopts.howtoperm                     = howtoperm;
regopts.tformrigid_allen_to_samp_20um = tformrigid;
regopts.howtoperm    = [3 1 2];
regopts.procpath     = sliceinfo.procpath;
regopts.registres    = sliceinfo.px_register;
regopts.allenres     = 10; % um
regopts.errall       = errall;
regopts.atlasaplims  = sliceinfo.atlasaplims;
regopts.pxsizes      = [sliceinfo.slicethickness/sliceinfo.px_register 1 1];
regopts.extentfactor = 6; % # slices to extend beyond rigid registration
save(fullfile(sliceinfo.procpath, 'regopts.mat'), '-struct', 'regopts')
%--------------------------------------------------------------------------
scalesize = [ceil(size(slicevol,[1 2])*sliceinfo.px_process/sliceinfo.px_register) Nslices];
volsave   = gather(squeeze(slicevol(:,:,1,:)));

if sliceinfo.medianfiltreg 
    for islice = 1:Nslices
        volsave(:, :, islice) = medfilt2(volsave(:, :, islice), medpx*[1 1]);
    end
end
volsave   = single(imresize3(volsave, scalesize));
regvolfac = (2^16-1)/max(volsave, [],[1 2]);
voldown   = uint16(regvolfac.*volsave);

samplepath = fullfile(sliceinfo.procpath, sprintf('sample_register_%dum.tif', regopts.registres));
options.compress = 'lzw';
options.message  = false;
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(voldown, samplepath, options);
%--------------------------------------------------------------------------
volsamp = single(permute(volsave, howtoperm));
volsamp = volsamp./quantile(volsamp, 0.999, 'all');
rout    = imref3d([Nslices size(avreg, [2 3])], 1, regopts.pxsizes(1), 1);
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
scalesize    = [ceil(size(slicevol,[1 2])*sliceinfo.px_process/sliceinfo.px_register) Nslices];
voldown      = zeros([scalesize(1:2) 3 scalesize(3)], 'uint8');

for ichan = 1:min(Nchans, 3)
    volproc     = imresize3(squeeze(slicevol(:, :, ichan, :)), scalesize);
    backproc    = single(median(sliceinfo.backvalues(ichan, :)));
    if backproc > 0
        volproc     = (single(volproc) - backproc)./backproc;
    else
        volproc     = single(volproc);
    end
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

