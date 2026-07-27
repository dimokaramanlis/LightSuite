function [slicevol, regopts] = alignSliceVolume(slicevol, sliceinfo)
%ALIGNSLICEVOLUME Align slice volume to Allen Atlas.
%   Processes one channel at a time to keep peak memory low:
%   channel 1 is loaded for point-cloud registration, then every channel
%   is individually loaded, reordered/centered, rigidly aligned and saved.

if isnumeric(slicevol)
    error('alignSliceVolume: numeric volume input is no longer supported. Pass the folder path (sliceinfo.slicevol).');
end
assert(isstring(slicevol) || ischar(slicevol), 'First argument must be a folder path string.');
slicevolpath = char(slicevol);

%--------------------------------------------------------------------------
% Read ordering/centering decisions — no image data needed at this point
fprintf('Reading ordering decisions... ');
orderfile      = fullfile(sliceinfo.procpath, 'volume_for_ordering_processing_decisions.txt');
tabledecisions = [];
sliceorder     = 1:sliceinfo.Nslices;
flipsdo        = false(sliceinfo.Nslices, 1);
toremove       = false(sliceinfo.Nslices, 1);
if exist(orderfile, 'file')
    tabledecisions = readtable(orderfile);
    sliceorder     = tabledecisions.NewOrderOriginalIndex;
    flipsdo        = tabledecisions.FlipState == 1;
    toremove       = tabledecisions.FlipState == -1;
end
toremove_ord = toremove(sliceorder);   % logical mask in post-reorder positions
fprintf('Done!\n');

%--------------------------------------------------------------------------
% Pre-compute centering parameters from the decisions table.
% Uses sliceinfo.size_proc (set by generateSliceVolume) — no image loading.
Ny         = sliceinfo.size_proc(1);
Nx         = sliceinfo.size_proc(2);
proc_scale = sliceinfo.px_register / sliceinfo.px_process;
cropsugg   = [];
pad_target = sliceinfo.size_proc;   % fallback when no centering suggestions exist

hasCentering = ~isempty(tabledecisions) && ...
    all(ismember({'xmin','xmax','ymin','ymax'}, tabledecisions.Properties.VariableNames));

if hasCentering
    cropsugg = round(double([tabledecisions.xmin, tabledecisions.xmax, ...
                              tabledecisions.ymin, tabledecisions.ymax]) * proc_scale);
    % mirror x limits for slices that were flipped left-right
    for iorig = find(flipsdo)'
        old_xmin = cropsugg(iorig, 1); old_xmax = cropsugg(iorig, 2);
        cropsugg(iorig, 1) = Nx - old_xmax + 1;
        cropsugg(iorig, 2) = Nx - old_xmin + 1;
    end
    cropsugg = cropsugg(sliceorder, :);
    cropsugg(toremove_ord, :) = [];
    cropsugg(:, [1 2]) = max(1, min(Nx, cropsugg(:, [1 2])));
    cropsugg(:, [3 4]) = max(1, min(Ny, cropsugg(:, [3 4])));
    crop_w     = cropsugg(:, 2) - cropsugg(:, 1) + 1;
    crop_h     = cropsugg(:, 4) - cropsugg(:, 3) + 1;
    pad_target = [max(crop_h), max(crop_w)];
    sliceinfo.size_proc = pad_target;
end

Nslices_now = sum(~toremove_ord);
Nchans      = numel(sliceinfo.channames);

%--------------------------------------------------------------------------
% Shared parameters
howtoperm = [3 1 2];
ireg      = 1;
finsize   = ceil(pad_target * sliceinfo.px_process / sliceinfo.px_register);
sizedown  = [finsize, Nslices_now];
medpx     = 2*floor((21./sliceinfo.px_process)/2) + 1;

%--------------------------------------------------------------------------
% Load and prepare REGISTRATION CHANNEL only (channel ireg).
% This channel stays in memory through the atlas registration step so it
% can be reused as the first iteration of the per-channel save loop.
fprintf('Loading registration channel... '); tic;
volchan1 = loadCenterChannel(slicevolpath, ireg, sliceorder, flipsdo, toremove_ord, cropsugg, pad_target);
fprintf('Done! Took %2.2f s\n', toc);

%--------------------------------------------------------------------------
fprintf('Standardizing and filtering volume... '); tic;
if sliceinfo.use_gpu
    volregister = single(gpuArray(imresize3(volchan1, sizedown)));
else
    volregister = single(imresize3(volchan1, sizedown));
end
fprintf('Done! Took %2.2f s\n', toc);

%==========================================================================
minperc = 0.25; maxperc = 0.75;
ny = sizedown(1); nx = sizedown(2);
xlook = round(nx*minperc):round(nx*maxperc);
ylook = round(ny*minperc):round(ny*maxperc);
allcoords = cell(Nslices_now, 1);

for islice = 1:Nslices_now
    currslice = volregister(:, :, islice);
    backval   = quantile(currslice(currslice>0), 0.01, 'all');
    dfimg     = (currslice - backval) ./ backval;
    dfimg(isinf(dfimg) | dfimg < 0) = 0;
    centvals  = dfimg(ylook, xlook);
    thresuse  = max(0.5, quantile(centvals, 0.99, 'all') / 4);
    ipx       = find(dfimg(:) > thresuse);
    [row, col] = ind2sub(size(dfimg), ipx);
    pccloud   = pointCloud(gather([col, row, zeros(size(row))]));
    Npts      = pccloud.Count;
    targetnum = min(2e4, Npts);
    pccloud   = pcdownsample(pccloud, 'random', targetnum/Npts);
    pccloud   = pcdenoise(pccloud, 'NumNeighbors', 100);
    allcoords{islice} = [pccloud.Location(:,1:2), islice * ones(pccloud.Count, 1)];
end
clear volregister;

allpts   = cat(1, allcoords{:});
yvals    = allpts(:,3) * sliceinfo.slicethickness / sliceinfo.px_register;
xvals    = allpts(:,2);
zvals    = allpts(:,1);
pcsample = pointCloud([xvals, yvals, zvals]);
%==========================================================================

fprintf('Loading and processing Allen Atlas template... '); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv        = niftiread(fullfile(allen_atlas_path, 'average_template_10.nii.gz'));
av        = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
tv        = tv(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
av        = av(sliceinfo.atlasaplims(1):sliceinfo.atlasaplims(2), :, :);
tvreg     = imresize3(tv, sliceinfo.px_atlas/sliceinfo.px_register);
avreg     = imresize3(av, sliceinfo.px_atlas/sliceinfo.px_register);
atlasframe = size(tv, [2 3]);
clear tv;

thresuse            = quantile(tvreg(tvreg>0), 0.99, 'all') / 4;
ipx                 = find(tvreg(:) > thresuse);
[row, col, lastdim] = ind2sub(size(tvreg), ipx);
tv_cloud = pointCloud(gather([col, row, lastdim]));
clear tvreg;
fprintf('Done! Took %2.2f s\n', toc);

%==========================================================================
tformslices(Nslices_now, 1) = rigidtform2d;
transformtypes = {'rigid', 'affine'};
transformsteps = [1 1 1];
errall = nan(numel(transformsteps), 1);
for istep = 1:numel(transformsteps)
    currtranstype               = transformtypes{transformsteps(istep)};
    fprintf('Optimization step %d/%d: %s\n', istep, numel(transformsteps), currtranstype);
    [tformrigid, errall(istep)] = alignAtlasToSample(tv_cloud, pcsample, tformslices);
    tformslices                 = refineSampleFromAtlas(tv_cloud, pcsample, tformrigid, currtranstype);
    [rrx, rry, rrz]             = reportRotationAngles(tformrigid.R);
    fprintf('%s\n', repmat('=', [1 75]));
end
%--------------------------------------------------------------------------
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
regopts.extentfactor = 10;
save(fullfile(sliceinfo.procpath, 'regopts.mat'), '-struct', 'regopts');

%--------------------------------------------------------------------------
% Pre-compute inspection-volume output size (matches original scalesize formula)
pxsamp         = sliceinfo.px_process / sliceinfo.px_register;
pxatlas        = sliceinfo.px_atlas   / sliceinfo.px_register;
aln_size       = ceil(atlasframe * pxatlas / pxsamp);      % spatial dims of aligned output
scalesize_insp = [ceil(aln_size * pxsamp), Nslices_now];   % downsampled to register res

%--------------------------------------------------------------------------
% Per-channel loop: load → reorder/center → align → save.
% Channel ireg (1) is already loaded; all others are loaded fresh.
fprintf('\nAligning and saving channels...\n');
makeNewDir(sliceinfo.slicevolfin);
opt_vol.compress = 'lzw'; opt_vol.message  = false;
opt_vol.color    = false;  opt_vol.big      = true;

dpsavelowres = fullfile(sliceinfo.procpath, 'volume_for_inspection.tiff');
voldown_insp = zeros([scalesize_insp(1:2), 3, scalesize_insp(3)], 'uint8');

for ci = 1:Nchans
    fprintf('Channel %d/%d: ', ci, Nchans);

    % --- load & center ---
    if ci == ireg
        fprintf('centering (pre-loaded)... ');
        volchan = volchan1;
        clear volchan1;
    else
        fprintf('loading and centering... '); tic;
        volchan = loadCenterChannel(slicevolpath, ci, sliceorder, flipsdo, toremove_ord, cropsugg, pad_target);
        fprintf('Done! Took %2.2f s  ', toc);
    end

    % --- rigid alignment ---
    fprintf('aligning... '); tic;
    sliceinfo_1ch            = sliceinfo;
    sliceinfo_1ch.backvalues = sliceinfo.backvalues(ci, :);
    volchan4d = reshape(volchan, [size(volchan,1), size(volchan,2), 1, Nslices_now]);
    clear volchan;
    finalchan = getRigidlyAlignedVolume(sliceinfo_1ch, volchan4d, tformslices, atlasframe);
    clear volchan4d;
    fprintf('Done! Took %2.2f s\n', toc);

    % --- registration TIFF + diagnostic plots for channel ireg ---
    if ci == ireg
        fprintf('  Saving registration TIFF... '); tic;
        scalesize_reg = [ceil(size(finalchan,[1 2]) * sliceinfo.px_process/sliceinfo.px_register), Nslices_now];
        volsave_reg   = gather(squeeze(finalchan));
        if sliceinfo.medianfiltreg
            for islice = 1:Nslices_now
                volsave_reg(:,:,islice) = medfilt2(volsave_reg(:,:,islice), medpx*[1 1]);
            end
        end
        volsave_reg = single(imresize3(volsave_reg, scalesize_reg));
        regvolfac   = (2^16-1) ./ max(volsave_reg, [], [1 2]);
        voldown_reg = uint16(regvolfac .* volsave_reg);
        samplepath  = fullfile(sliceinfo.procpath, sprintf('sample_register_%dum.tif', regopts.registres));
        opt_tif.compress = 'lzw'; opt_tif.message = false;
        if exist(samplepath, 'file'), delete(samplepath); end
        saveastiff(voldown_reg, samplepath, opt_tif);
        clear voldown_reg;
        fprintf('Done! Took %2.2f s\n', toc);

        fprintf('  Generating diagnostic plots... ');
        volsamp = single(permute(volsave_reg, howtoperm));
        volsamp = volsamp ./ quantile(volsamp, 0.999, 'all');
        clear volsave_reg;
        rout   = imref3d([Nslices_now, size(avreg, [2 3])], 1, regopts.pxsizes(1), 1);
        avtest = imwarp(avreg, imref3d(size(avreg)), tformrigid, 'nearest', 'OutputView', rout);
        for idim = 1:3
            cf = plotAnnotationComparison(uint8(255*volsamp), avtest, idim, regopts.pxsizes);
            print(cf, fullfile(sliceinfo.procpath, ...
                sprintf('%s_dim%d_initial_registration', sliceinfo.mousename, idim)), '-dpng');
            close(cf);
        end
        clear volsamp;
        fprintf('Done!\n');
    end

    % --- save aligned channel to slicevolfin ---
    fprintf('  Saving channel %d... ', ci); tic;
    currname = sprintf('chan%02d_%s.tiff', ci, sliceinfo.channames{ci});
    volpath  = fullfile(sliceinfo.slicevolfin, currname);
    if exist(volpath, 'file'), delete(volpath); end
    saveastiff(squeeze(finalchan), volpath, opt_vol);
    fprintf('Done! Took %2.2f s\n', toc);

    % --- contribute to inspection volume (first 3 channels) ---
    if ci <= 3
        volproc  = single(imresize3(squeeze(finalchan), scalesize_insp));
        clear finalchan;
        backproc = single(median(sliceinfo.backvalues(ci, :)));
        if backproc > 0
            volproc = (volproc - backproc) ./ backproc;
        end
        maxval = quantile(volproc, 0.999, 'all');
        voldown_insp(:,:,ci,:) = uint8(255 * volproc / maxval);
        clear volproc;
    else
        clear finalchan;
    end
end

%--------------------------------------------------------------------------
fprintf('Saving inspection volume... '); tic;
opt_insp.compress = 'lzw'; opt_insp.message = false;
opt_insp.color    = true;  opt_insp.big     = false;
if exist(dpsavelowres, 'file'), delete(dpsavelowres); end
saveastiff(voldown_insp, dpsavelowres, opt_insp);
clear voldown_insp;
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------

slicevol = [];  % volume is saved to disk; not returned as array
end

% =========================================================================
function volchan = loadCenterChannel(slicevolpath, ci, sliceorder, flipsdo, toremove_ord, cropsugg, pad_target)
%LOADCENTERCHANNEL Load one channel, apply reorder/flip/remove, apply centering.
%   Returns [Ny, Nx, Nslices_now] uint16.

volchan = loadLargeSliceVolume(slicevolpath, ci);   % squeeze inside: [Ny, Nx, Nslices]
if ndims(volchan) == 4                              % guard: single-chan may still be 4-D
    volchan = squeeze(volchan(:, :, 1, :));
end

% reorder / flip / remove (same sequence as original code)
volchan(:, :, flipsdo)      = flip(volchan(:, :, flipsdo), 2);
volchan                     = volchan(:, :, sliceorder);
volchan(:, :, toremove_ord) = [];

if isempty(cropsugg); return; end   % no centering requested

% centering: crop to brain region, then pad to uniform pad_target
Nslices_now  = size(volchan, 3);
slicevol_cen = zeros([pad_target, Nslices_now], 'uint16');
for islice = 1:Nslices_now
    xmin_s = cropsugg(islice, 1); xmax_s = cropsugg(islice, 2);
    ymin_s = cropsugg(islice, 3); ymax_s = cropsugg(islice, 4);
    padsz   = pad_target - [ymax_s - ymin_s + 1, xmax_s - xmin_s + 1];
    padleft = floor(padsz / 2);
    patch   = volchan(ymin_s:ymax_s, xmin_s:xmax_s, islice);
    backval = uint16(quantile(single(patch(patch > 0)), 0.01, 'all'));
    patch   = padarray(patch, padleft, backval, 'pre');
    patch   = padarray(patch, padsz - padleft, backval, 'post');
    slicevol_cen(:, :, islice) = patch;
end
volchan = slicevol_cen;
end
