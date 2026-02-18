function [slicenew, regopts] = coarseAlignSliceVolumeToAtlas(sliceinfo, slicevol, fillvalues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rng(1);
regopts.howtoperm = [3 1 2];
regopts.procpath  = sliceinfo.procpath;
regopts.registres = sliceinfo.px_register;
regopts.allenres  = 10; % um
[Nchan, Nslices]  = size(slicevol, [3 4]);
pxthick           = sliceinfo.slicethickness/sliceinfo.px_register;
%--------------------------------------------------------------------------
fprintf('Loading and processing Allen Atlas template... '); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tvreg            = imresize3(tv, regopts.allenres/sliceinfo.px_register);
avreg            = imresize3(av, regopts.allenres/sliceinfo.px_register, "Method", "nearest", 'AntiAliasing',true);
limskeep         = [55, size(tvreg, 1)-100]; % exclude cerebellum and olfactory
tv_cloud         = extractHighSFVolumePoints(tvreg, sliceinfo.px_register, limskeep);
Npts             = tv_cloud.Count;
tv_cloud_use     = pcdownsample(tv_cloud,'random', 10000/Npts, 'PreserveStructure',true);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Extracting points of interest from data... '); tic;

alignedvol  = squeeze(slicevol(:, :, 1, :));
scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
alignedvol  = single(alignedvol);
bval        = median(single(sliceinfo.backvalues(1,:)));
alignedvol  = (alignedvol - bval)/bval;
alignedvol(alignedvol<0) = 0;
alignedvol = imresize3(alignedvol, [finsize Nslices]);


% we perform spatial bandpass filtering
imhigh          = spatial_bandpass(alignedvol, scalefilter, 3, 3, sliceinfo.use_gpu);
thresuse        = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp           = find(imhigh>thresuse);
[row,col,slice] = ind2sub(size(imhigh), idxcp);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% we first align 3d volume to Allen
fprintf('Coarse alignment of slice volume to Allen atlas... '); tic;
zvals = slice * pxthick;
xvals = col ;
yvals = row;
Xinput = [gather(yvals), gather(zvals), gather(xvals)];
% Xinput = Xinput - mean(Xinput);
pcall  = pointCloud(Xinput);
pcplot = pcdownsample(pcall, 'random', 0.05, 'PreserveStructure',true);

[tformrigid, pcreginit, resinit] = pcregistercpd(tv_cloud_use, pcplot, "Transform","Rigid",...
    "Verbose",false,"OutlierRatio",0.00, 'MaxIterations', 150, 'Tolerance', 1e-6);
tvtrans = pctransform(tv_cloud, tformrigid);
% figure;
% pcshowpair(pcplot, pcreg)

fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% based on the first alignment, we transform each slice so that it matches
% the corresponding Allen slice in 2d
%%
facuse        = 1.1; % compensation factor for elongated slices that may end up outside atlas boundaries
sizeout       = round(facuse*ceil(size(tvreg,[2 3])*sliceinfo.px_register/sliceinfo.px_process));
R_out_fullres = imref2d(sizeout);
Rsample       = imref2d(size(slicevol,[1 2]));
slicenew      = zeros([R_out_fullres.ImageSize Nchan Nslices], 'uint16');
sigmause      = pxthick/4;

for islicecurr = 1:sliceinfo.Nslices

    fprintf('Aligning slice %02d to atlas... ', islicecurr); 
    slicetic = tic;

    % we find the x-y points in the original atlas space
    Xslicecurr = Xinput(slice == islicecurr, :);
    apcurr     = Xslicecurr(1,2);
    iatlascurr = abs(tvtrans.Location(:,2) - apcurr) < sliceinfo.slicethickness/sliceinfo.px_register/2;
    Xatlascurr = tv_cloud.Location(iatlascurr, :);

    Xslicecurr(:, 2) = randn(size(Xslicecurr, 1), 1) *sigmause;
    Xatlascurr(:, 2) = randn(size(Xatlascurr, 1), 1) *sigmause;

    pcatlascurr = pointCloud(facuse*Xatlascurr*sliceinfo.px_register/sliceinfo.px_process);
    pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample',50, 'PreserveStructure',true);

    pcslicecurr = pointCloud(Xslicecurr*sliceinfo.px_register/sliceinfo.px_process);
    pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', 6, 'PreserveStructure',true);

    [R,T,data2] = icp(pcslicecurr.Location(:,[3 1]), pcatlascurr.Location(:,[3 1]), 100, 10, 1, 1e-6);
    res         = mean(min(pdist2(data2', pcslicecurr.Location(:,[3 1])),[],1))*sliceinfo.px_process;

    tformcurr   = rigidtform2d(R, T);
    tformf      = tformcurr.invert;


    for ichan = 1:Nchan
        slicenew(:,:,ichan, islicecurr) =  imwarp(slicevol(:, :, ichan, islicecurr), Rsample, tformf, ...
                'linear', 'OutputView', R_out_fullres, 'FillValues', fillvalues(ichan));
    end

    fprintf('Done! Error %2.2f um. Took %2.2f s\n', res, toc(slicetic));
    %-----------------------------------------------------------------------
    % datatest = tformcurr.transformPointsInverse(pcslicecurr.Location(:,[3 1]));
    % [tformcod, regcpd] = pcregistercpd(pcslicecurr, pcatlascurr, ...
    %     'Transform','Rigid', 'Tolerance',1e-6,'OutlierRatio',0.0,'Verbose',true,'MaxIterations',100);
    % 
    % clf;
    % subplot(1,3,1)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % subplot(1,3,2)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     datatest(:,1), datatest(:,2), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % title('icp')
    % subplot(1,3,3)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     regcpd.Location(:,3), regcpd.Location(:,1), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % title('cpd')
    % pause;
    %-----------------------------------------------------------------------
end
%%
%--------------------------------------------------------------------------
% we then re-align in 3d the proper volume to Allen

fprintf('Initializing atlas registration...\n'); tic;

alignedvol  = squeeze(slicenew(:, :, 1, :));
scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(size(alignedvol,[1 2])*sliceinfo.px_process/sliceinfo.px_register);
alignedvol  = single(alignedvol);
bval        = median(single(sliceinfo.backvalues(1,:)));
alignedvol  = (alignedvol - bval)/bval;
alignedvol(alignedvol<0) = 0;
alignedvol = imresize3(alignedvol, [finsize Nslices]);

% we perform spatial bandpass filtering
imhigh          = spatial_bandpass(alignedvol, scalefilter, 3, 3, sliceinfo.use_gpu);
thresuse        = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp           = find(imhigh>thresuse);
[row,col,slice] = ind2sub(size(imhigh), idxcp);

zvals = (slice * sliceinfo.slicethickness)/sliceinfo.px_register;
xvals = col ;
yvals = row;
Xinput = [gather(yvals), gather(zvals), gather(xvals)];
% Xinput = Xinput - mean(Xinput);
pcall  = pointCloud(Xinput);
pcplot = pcdownsample(pcall, 'random', 0.05, 'PreserveStructure',true);

[tformrigidfin, pcregfin, resfin] = pcregistercpd(tv_cloud_use, pcplot, "Transform","Rigid",...
    "Verbose",false,"OutlierRatio",0.00, 'MaxIterations', 150, 'Tolerance', 1e-6);
fprintf('Original error: %2.3f. Final error: %2.3f.\n', resinit, resfin); 
%--------------------------------------------------------------------------
% save volume for control point and registration

regvolfac = (2^16-1)/max(alignedvol, [],"all");
voldown   = uint16(regvolfac*alignedvol);

samplepath = fullfile(sliceinfo.procpath, sprintf('sample_register_%dum.tif', regopts.registres));
options.compress = 'lzw';
options.message  = false;
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(voldown, samplepath, options);
%--------------------------------------------------------------------------
% let's plot results

volsamp  = permute(alignedvol, regopts.howtoperm);
volsamp  = volsamp./quantile(volsamp, 0.999, "all");
raatlas  = imref3d(size(avreg),  1, 1, 1);
rasample = imref3d(size(volsamp), 1, sliceinfo.slicethickness/ sliceinfo.px_register, 1);
avtest   = imwarp(avreg, raatlas, tformrigidfin,'nearest', 'OutputView',rasample);

pxsizes = [rasample.PixelExtentInWorldY rasample.PixelExtentInWorldX rasample.PixelExtentInWorldZ];
for idim = 1:3
    cf = plotAnnotationComparison(uint8(255*volsamp), avtest, idim, pxsizes);
    print(cf, fullfile( sliceinfo.procpath, sprintf('dim%d_initial_registration', idim)), '-dpng')
    close(cf);
end

regopts.tform_rigid_AllenToSample_20um = tformrigidfin;
regopts.pxsizes                        = pxsizes;
save(fullfile(sliceinfo.procpath, 'regopts.mat'), '-struct', 'regopts')
%--------------------------------------------------------------------------
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
end