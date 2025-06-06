function [slicenew, tformAllenToSlice] = coarseAlignSliceVolumeToAtlas(sliceinfo, slicevol, fillvalues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Nchan, Nslices] = size(slicevol, [3 4]);
%--------------------------------------------------------------------------
fprintf('Loading and processing Allen Atlas template... '); tic;
allenres         = 10; % um
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
tvreg            = imresize3(tv, allenres/sliceinfo.px_register);
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
zvals = (slice * sliceinfo.slicethickness)/sliceinfo.px_register;
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

R_out_fullres = imref2d(ceil(size(tvreg,[2 3])*sliceinfo.px_register/sliceinfo.px_process));
Rsample       = imref2d(size(slicevol,[1 2]));
slicenew      = zeros([R_out_fullres.ImageSize Nchan Nslices], 'uint16');


for islicecurr = 1:sliceinfo.Nslices

    fprintf('Aligning slice %02d to atlas... ', islicecurr); 
    slicetic = tic;

    % we find the x-y points in the original atlas space
    Xslicecurr = Xinput(slice == islicecurr, :);
    apcurr     = Xslicecurr(1,2);
    iatlascurr = abs(tvtrans.Location(:,2) - apcurr) < sliceinfo.slicethickness/sliceinfo.px_register/2;


    pcatlascurr = pointCloud(tv_cloud.Location(iatlascurr, :)*sliceinfo.px_register/sliceinfo.px_process);
    pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample',50, 'PreserveStructure',true);

    pcslicecurr = pointCloud(Xslicecurr*sliceinfo.px_register/sliceinfo.px_process);
    pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', 6, 'PreserveStructure',true);

    [R,T,data2] = icp(pcslicecurr.Location(:,[3 1]), pcatlascurr.Location(:,[3 1]), 100, 10, 1, 1e-6);
    res         = mean(min(pdist2(data2', pcslicecurr.Location(:,[3 1])),[],1))*sliceinfo.px_process;

    tformcurr   = rigidtform2d(R, T);
    tformf      = tformcurr.invert;
    

    % testim = imwarp(slicevol(:, :, 1, islicecurr), Rsample, tformf, ...
    %         'linear', 'OutputView', R_out_fullres, 'FillValues', 0);

    for ichan = 1:Nchan
        slicenew(:,:,ichan, islicecurr) =  imwarp(slicevol(:, :, ichan, islicecurr), Rsample, tformf, ...
                'linear', 'OutputView', R_out_fullres, 'FillValues', fillvalues(ichan));
    end

    % pcatlascurr = pointCloud(tvtrans.Location(iatlascurr, :));
    % pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample',50, 'PreserveStructure',true);
    % pcslicecurr = pointCloud(Xslicecurr);
    % pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', 6, 'PreserveStructure',true);
    % 
    % [R,T,data2] = icp(pcslicecurr.Location(:,[1 3]), pcatlascurr.Location(:,[1 3]), 100, 10, 1, 1e-6);
    % tformcurr = rigidtform2d(R', flip(T));
    % tformf    = tformcurr.invert;
    % tformsave(islicecurr) = rigidtform2d(tformf.A);
    % data2 = data2';
    % res = mean(min(pdist2(data2, pcslicecurr.Location(:,[1 3])),[],1));

    fprintf('Done! Error %2.2f um. Took %2.2f s\n', res, toc(slicetic));
    
    % datatest = tformcurr.transformPointsInverse(pcslicecurr.Location(:,[1 3]));


    % clf;
    % subplot(2,2,1)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % subplot(2,2,2)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     datatest(:,2), datatest(:,1), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % 
    % plot(data2(:,2), data2(:,1), 'r.',...
    %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
    % newvol(:,:,islicecurr) = imwarp(alignedvol(:,:,islicecurr), rasample_2d_lowres,  tformsave(islicecurr), ...
    %     'linear', 'OutputView', raatlas_2d_lowres);
    % 
    % T_rigid = tformsave(islicecurr);
    % T_composite = affinetform2d(T_rigid.A*T_scale.A );

    % subplot(2,2,3)
    % % imagesc(alignedvol(:,:,islicecurr));hold on;
    % % plot(pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'r.')
    % % axis equal; axis tight;
    % imshowpair(squeeze(tvinsamp(islicecurr, :, :)), alignedvol(:,:,islicecurr))
    % subplot(2,2,4)
    % % imagesc(resamp(:,:,islicecurr));hold on;
    % % plot(pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'r.')
    % % axis equal; axis tight;
    % imshowpair(squeeze(tvinsamp(islicecurr, :, :)),resamp)
    % pause;
end
%%
%--------------------------------------------------------------------------
% we then re-align in 3d the proper volume to Allen
% 
% fprintf('Extracting points of interest from data... '); tic;
% 
% alignedvol  = squeeze(slicenew(:, :, 1, :));
% scalefilter = 100/sliceinfo.px_register;
% finsize     = ceil(size(alignedvol,[1 2])*sliceinfo.px_process/sliceinfo.px_register);
% alignedvol  = single(alignedvol);
% bval        = median(single(sliceinfo.backvalues(1,:)));
% alignedvol  = (alignedvol - bval)/bval;
% alignedvol(alignedvol<0) = 0;
% alignedvol = imresize3(alignedvol, [finsize Nslices]);
% 
% 
% % we perform spatial bandpass filtering
% imhigh          = spatial_bandpass(alignedvol, scalefilter, 3, 3, sliceinfo.use_gpu);
% thresuse        = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
% idxcp           = find(imhigh>thresuse);
% [row,col,slice] = ind2sub(size(imhigh), idxcp);
% fprintf('Done! Took %2.2f s\n', toc); 
% 
% 
% zvals = (slice * sliceinfo.slicethickness)/sliceinfo.px_register;
% xvals = col ;
% yvals = row;
% Xinput = [gather(yvals), gather(zvals), gather(xvals)];
% % Xinput = Xinput - mean(Xinput);
% pcall  = pointCloud(Xinput);
% pcplot = pcdownsample(pcall, 'random', 0.05, 'PreserveStructure',true);
% 
% [tformrigidfin, pcregfin, resfin] = pcregistercpd(tv_cloud_use, pcplot, "Transform","Rigid",...
%     "Verbose",true,"OutlierRatio",0.00, 'MaxIterations', 150, 'Tolerance', 1e-6);
tformAllenToSlice = tformrigid;
%--------------------------------------------------------------------------
end