function opts = initializeRegistration(opts, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function re-orients the sample to match the atlas and facilitate
% control point selection
%==========================================================================
if nargin < 2
    backvol  = readDownStack(opts.regvolpath);
else
    backvol = varargin{1};
end
downfac = opts.atlasres/opts.registres;
%==========================================================================
%%
% we fist preprocess the registration volume
centpx    = round(size(backvol)/2);

naround   = round(min(centpx)/3);
centind   = round(size(backvol)/2) + (-naround:naround)';
bvalinds  = randperm(numel(backvol), 1e5);
topval    = quantile(backvol(centind(:,1), centind(:,2),centind(:,3)), 0.999,'all')*2;
% bottomval = quantile(backvol(bvalinds), 0.01, 'all');
bottomval = 0;
newvol    = (single(backvol)-single(bottomval))/single(topval-bottomval);

% newvol    = uint8(newvol * 255);
%==========================================================================
% we then prepare the volume and extract corresponding points 
% volumereg  = imresize3(newvol, 0.5);
volumereg  = permute(newvol, [1 3 2]);
ls_cloud   = extractVolumePoints(uint16(volumereg * (2^16-1)), 15);

%==========================================================================
% load atlas and extract corresponding points
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
% atlas needs to match volume dimensions
% tv    = permute(tv, [1 3 2]);
tvreg = imresize3(tv, downfac);
avreg = imresize3(av, downfac, 'Method','nearest');

tv_cloud = extractVolumePoints(tvreg, 15);

% tvaffine = imwarp(tv, Rmoving, tform, 'OutputView',Rfixed);
%==========================================================================
% the first step is to make sure our sample is nicely aligned

Rtemplate = imref3d(size(tvreg));
Rsample   = imref3d(size(volumereg));

transinit = pcregistercpd(ls_cloud, tv_cloud,'Transform', 'Rigid','OutlierRatio',0.0,...
    'Verbose',false,'MaxIterations',100,'Tolerance',1e-6);
% volumetest = imwarp(volumereg, Rsample, transinit, 'OutputView',Rtemplate);
%==========================================================================
% we plot the first step
% viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
% volshow(volumetest,Parent=viewerUnregistered,RenderingStyle="Isosurface", Colormap=[1 0 1],Alphamap=1);
% volshow(tvreg,Parent=viewerUnregistered,RenderingStyle="Isosurface", Colormap=[0 1 0],Alphamap=1);

avtest = imwarp(avreg, Rtemplate, transinit.invert,'nearest', 'OutputView',Rsample);
volmax = quantile(volumereg,0.999,'all');
for idim = 1:3
    cf = plotAnnotationComparison(uint8(255*volumereg/volmax), avtest, idim);
    savepngFast(cf, opts.savepath, sprintf('dim%d_initial_registration', idim), 300, 2);
    close(cf);
end
%==========================================================================
% let's save stuff
opts.permute_sample_to_atlas = [1 3 2];
opts.original_trans          = transinit;
opts.downfac_reg             = downfac;

% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')
%==========================================================================
end