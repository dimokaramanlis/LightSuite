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

volregister = imresize3(volregister, sizedown);

if sliceinfo.use_gpu
    volregister = single(gpuArray(volregister));
else
    volregister = single(volregister);
end
fprintf('Done! Took %2.2f s\n', toc); 
%==========================================================================
minperc   = 0.25;
maxperc   = 0.75;
ny = sizedown(1);
nx = sizedown(2);
xlook = round(nx*minperc):round(nx*maxperc);
ylook = round(ny*minperc):round(ny*maxperc);

allcoords = cell(Nslices, 1);

for islice = 1:Nslices
    currslice = volregister(:, :, islice);
    backval = quantile(currslice(currslice>0), 0.01, 'all');
    dfimg  = (currslice- backval)./backval;
    dfimg(isinf(dfimg) | dfimg < 0) = 0;

    centvals = dfimg(ylook, xlook);
    thresuse = max(0.5, quantile(centvals,0.99,'all')/4);

    ipx = find(dfimg(:) > thresuse);
    [row, col] = ind2sub(size(dfimg), ipx);
    pccloud = pointCloud(gather([col, row, zeros(size(row))]));
    Npts    = pccloud.Count;
    targetnum = min(2e4, Npts);
    pccloud = pcdownsample(pccloud, 'random', targetnum/Npts);
    pccloud = pcdenoise(pccloud, 'NumNeighbors',100);
    allcoords{islice} = [pccloud.Location(:,1:2)  islice * ones(pccloud.Count, 1)];
end
allpts = cat(1, allcoords{:});
yvals      = allpts(:,3) * sliceinfo.slicethickness/sliceinfo.px_register;
xvals      = allpts(:,2);
zvals      = allpts(:,1);
Xinput     = [xvals, yvals, zvals];
pcsample   = pointCloud(Xinput);
%==========================================================================

fprintf('Loading and processing Allen Atlas template... '); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tv               = tv(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
av               = av(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
tvreg            = imresize3(tv, sliceinfo.px_atlas/sliceinfo.px_register);
avreg            = imresize3(av, sliceinfo.px_atlas/sliceinfo.px_register);
atlasframe = size(tv, [2 3]);

thresuse         = quantile(tvreg(tvreg>0),0.99,'all')/4;
ipx = find(tvreg(:) > thresuse);
[row, col, lastdim] = ind2sub(size(tvreg), ipx);
tv_cloud = pointCloud(gather([col, row, lastdim]));


% Npts             = size(tv_points, 1);
% tv_cloud_use     = pcdownsample(pointCloud(tv_points),'random', 20000/Npts, 'PreserveStructure',true);
fprintf('Done! Took %2.2f s\n', toc); 
%==========================================================================
%--------------------------------------------------------------------------
% optimization steps
tformslices(Nslices, 1) = rigidtform2d;


% Tvec       = median(pcsample.Location)- range(tv_cloud.Location)/2;
% tformrigid = rigidtform3d(eye(3),Tvec);
% % 
% tvclouddown = pcdownsample(tv_cloud, 'random', (2e4/tv_cloud.Count));
% pcclouddown = pcdownsample(pcsample, 'random', (2e4/pcsample.Count));
% 
% [tformrigid, movreg] = pcregisterndt(tvclouddown,pcclouddown,50, 'Verbose',true);

transformtypes = {'rigid', 'affine'};
transformsteps = [1 1 1];
errall = nan(numel(transformsteps), 1);
for istep = 1:numel(transformsteps)
    currtranstype               = transformtypes{transformsteps(istep)};
    fprintf('Optimization step %d/%d: %s\n', istep, numel(transformsteps), currtranstype)
    [tformrigid, errall(istep)] = alignAtlasToSample(tv_cloud, pcsample, tformslices);
    tformslices                 = refineSampleFromAtlas(tv_cloud, pcsample, tformrigid, currtranstype);
    [rrx, rry, rrz]             = reportRotationAngles(tformrigid.R);
    fprintf('%s\n', repmat('=', [1 75]));
end
%--------------------------------------------------------------------------

slicevol = getRigidlyAlignedVolume(sliceinfo, slicevol, tformslices, atlasframe);

regopts = struct();
regopts.howtoperm                     = howtoperm;
regopts.tformrigid_allen_to_samp_20um = tformrigid;
regopts.howtoperm    = [3 1 2];
regopts.procpath     = sliceinfo.procpath;
regopts.registres    = sliceinfo.px_register;
regopts.processres   = sliceinfo.px_process;
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
regvolfac = (2^16-1)./max(volsave, [],[1 2]);
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

