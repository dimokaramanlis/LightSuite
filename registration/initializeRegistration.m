function opts = initializeRegistration(opts, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function re-orients the sample to match the atlas and facilitate
% control point selection
%==========================================================================
if nargin < 2
    fidreg  = fopen(opts.regvolpath, 'r');
    backvol = fread(fidreg,[prod(opts.regvolsize) 1],"*uint16");
    fclose(fidreg);
    backvol = reshape(backvol, opts.regvolsize);
else
    backvol = varargin{1};
end

%==========================================================================
%%
% we fist preprocess the registration volume
centpx    = round(size(backvol)/2);

naround   = round(min(centpx)/3);
centind   = round(size(backvol)/2) + (-naround:naround)';
bvalinds  = randperm(numel(backvol), 1e5);
topval    = quantile(backvol(centind(:,1), centind(:,2),centind(:,3)), 0.999,'all')*2;
bottomval = quantile(backvol(bvalinds), 0.01, 'all');
newvol    = (single(backvol)-single(bottomval))/single(topval-bottomval);

% newvol    = uint8(newvol * 255);
%==========================================================================
% we then prepare the volume and extract corresponding points 
volumereg  = imresize3(newvol, 0.5);
volumereg  = permute(volumereg, [1 3 2]);

ls_cloud = extractVolumePoints(uint16(volumereg * (2^16-1)), 15);

%==========================================================================
% load atlas and extract corresponding points
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um.npy'));
% atlas needs to match volume dimensions
% tv    = permute(tv, [1 3 2]);
tvreg = imresize3(tv, 0.5);
avreg = imresize3(av, 0.5, 'Method','nearest');

tv_cloud = extractVolumePoints(tvreg, 15);

% tvaffine = imwarp(tv, Rmoving, tform, 'OutputView',Rfixed);
%==========================================================================
% the first step is to make sure our sample is nicely aligned

Rtemplate = imref3d(size(tvreg));
Rsample   = imref3d(size(volumereg));

transinit = pcregistercpd(ls_cloud, tv_cloud,'Transform', 'Rigid','OutlierRatio',0.0,...
    'Verbose',false,'MaxIterations',100,'Tolerance',1e-6);
volumetest = imwarp(volumereg, Rsample, transinit, 'OutputView',Rtemplate);
%==========================================================================
% we plot the first step
% viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
% volshow(volumetest,Parent=viewerUnregistered,RenderingStyle="Isosurface", Colormap=[1 0 1],Alphamap=1);
% volshow(tvreg,Parent=viewerUnregistered,RenderingStyle="Isosurface", Colormap=[0 1 0],Alphamap=1);

for idim = 1:3
    cf = plotAnnotationComparison(uint8(255*volumetest/0.5), single(avreg), idim);
    savepngFast(cf, opts.savepath, sprintf('dim%d_initial_registration', idim), 300, 2);
    close(cf);
end
%==========================================================================
% let's save stuff

samplepath = fullfile(opts.savepath, 'sample_register_20um.tif');

opts.permute_sample_to_atlas = [1 3 2];
opts.original_trans          = transinit;
opts.downfac_reg             = 0.5;
opts.samplepath_reg          = samplepath;
opts.regvolsize_down         = size(volumereg);

% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')

% save volumes for control points
options.compress = 'no';
volsave = uint8(volumetest*255);
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(volsave, samplepath, options);
%==========================================================================
end