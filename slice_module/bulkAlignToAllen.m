function regopts = bulkAlignToAllen(sliceinfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
regopts.howtoperm = [3 1 2];
%--------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
alignedvol = loadLargeSliceVolume(sliceinfo.slicevolfin, 1);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Loading and processing Allen Atlas template... '); tic;
allenres = 10;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
tvreg            = imresize3(tv, allenres/sliceinfo.px_register);
limskeep         = [55, size(tvreg, 1)-100]; % exclude cerebellum and olfactory
tv_cloud         = extractHighSFVolumePoints(tvreg, sliceinfo.px_register, limskeep);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Extracting points of interest from data... '); tic;

scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
alignedvol  = single(alignedvol);
bval        = median(single(sliceinfo.backvalues(1,:)));
alignedvol  = (alignedvol - bval)/bval;
alignedvol(alignedvol<0) = 0;
alignedvol = imresize3(alignedvol, [finsize sliceinfo.Nslices]);


% we perform spatial bandpass filtering
imhigh          = spatial_bandpass(alignedvol, scalefilter, 3, 3, sliceinfo.use_gpu);
thresuse        = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp           = find(imhigh>thresuse);
[row,col,slice] = ind2sub(size(imhigh), idxcp);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------

%%
zvals = (slice * sliceinfo.slicethickness)/sliceinfo.px_register;
xvals = col ;
yvals = row;
Xinput = [gather(yvals), gather(zvals), gather(xvals)];
% Xinput = Xinput - mean(Xinput);
pcall  = pointCloud(Xinput);
pcplot = pcdownsample(pcall, 'random', 0.05, 'PreserveStructure',true);
clf;
scatter3(pcplot.Location(:,1),pcplot.Location(:,2),pcplot.Location(:,3), 2, 'filled');
hold on;
scatter3(tv_cloud.Location(:,1),tv_cloud.Location(:,2),tv_cloud.Location(:,3), 2, 'filled')

%%
[tformrigid, pcreg] = pcregistercpd(tv_cloud, pcplot, "Transform","Rigid",...
    "Verbose",true,"OutlierRatio",0.00, 'MaxIterations', 100, 'Tolerance', 1e-6);
figure;
pcshowpair(pcplot, pcreg)
%%
volsamp  = permute(alignedvol, regopts.howtoperm);
raatlas  = imref3d(size(tvreg),  1, 1, 1);
rasample = imref3d(size(volsamp), 1, sliceinfo.slicethickness/ sliceinfo.px_register, 1);
tvinsamp = imwarp(tvreg, raatlas, tformrigid, 'OutputView', rasample);
%%
figure;
for ii = 1:sliceinfo.Nslices
    subplot(1,2,1); 
    imagesc(squeeze(tvinsamp(ii,:,:))); 
    axis equal; axis tight;
    subplot(1,2,2); 
    imagesc(squeeze(volsamp(ii,:,:))); 
    axis equal; axis tight;
    pause; 
end

end