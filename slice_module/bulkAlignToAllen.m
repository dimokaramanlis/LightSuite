function [outputArg1,outputArg2] = bulkAlignToAllen(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



allenres = 10;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
tvreg            = imresize3(tv, allenres/sliceinfo.px_register);
tvreg([1:50, end-80:end], :, :) = 0; % exclude cerebellum and olfactory
tv_cloud         = extractHighSFVolumePoints(tvreg, sliceinfo.px_register);


%%
ireg = 1;
scalefilter = 100/sliceinfo.px_register;
volregister = squeeze(alignedvol(:, :, ireg, :));
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
volregister = single(volregister);
bval        = median(single(sliceinfo.backvalues(ireg,:)));
volregister = (volregister - bval)/bval;
volregister(volregister<0) = 0;
volregister = imresize3(volregister, [finsize sliceinfo.Nslices]);

% we perform spatial bandpass filtering
imhigh         = spatial_bandpass(volregister, scalefilter, 3, 3, sliceinfo.use_gpu);
thresuse       = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp           = find(imhigh>thresuse);
[row,col,slice] = ind2sub(size(imhigh), idxcp);
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
pcshowpair(pcplot, pcreg)
end