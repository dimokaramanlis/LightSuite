function opts = initializeRegistration(inputpath, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function re-orients the sample to match the atlas and facilitate
% control point selection
%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) isstring(x) || ischar(x));
addParameter(p, 'FlipX', false, @islogical);
addParameter(p, 'Volume', [], @isnumeric);
addParameter(p, 'bpcdpath', 'C:\GitHub\bcpd\win\bcpd.exe', ...
    @(x) isstring(x) | ischar(x) );
parse(p, inputpath, varargin{:});
params = p.Results;
%==========================================================================
optsfile   = dir(fullfile(inputpath, '*regopts.mat'));
if numel(optsfile) > 1
    ikeep      = find(contains({optsfile.name}', '561'));
else
    ikeep      = 1;
end
if ~isempty(optsfile)
    opts = load(fullfile(optsfile(ikeep).folder, optsfile(ikeep).name));
    opts = opts.opts;
end
%==========================================================================
if ~isempty(params.Volume)
    backvol = params.Volume;
else
    backvol = readDownStack(opts.regvolpath);
end
downfac = opts.atlasres/opts.registres;
%==========================================================================
flipvec = [params.FlipX, false, false];
Tflip   = affinetform3d(createFlipTransform(size(backvol), flipvec));
%==========================================================================
%%
% we fist preprocess the registration volume
centpx    = round(size(backvol)/2);

naround   = round(min(centpx)/3);
centind   = round(size(backvol)/2) + (-naround:naround)';
topval    = quantile(backvol(centind(:,1), centind(:,2),centind(:,3)), 0.999,'all')*2;
bottomval = 0;
newvol    = (single(backvol)-single(bottomval))/single(topval-bottomval);
%==========================================================================
% we then prepare the volume and extract corresponding points 
% volumereg  = imresize3(newvol, 0.5);
fprintf('Creating cloud for sample volume... '); tic;
volumereg  = permute(newvol, [1 3 2]);
ls_cloud   = extractSamplePoints(volumereg, 5);
fprintf('Done! Took %2.1f s. Found %d points.\n', toc, ls_cloud.Count);
%==========================================================================
% load atlas and extract corresponding points
fprintf('Loading atlas data and generating the atlas cloud... '); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tvreg   = imresize3(tv, downfac);
avreg   = imresize3(av, downfac, 'Method','nearest');

% tv_cloud = extractVolumePoints(tvreg, 15);
tvforpoints             = single(tvreg);
tvforpoints(avreg == 0) = 0;
tv_cloud                = extractVolumePointsGradient(tvforpoints, 20, 5);
fprintf('Done! Took %2.1f s. Found %d points.\n', toc, tv_cloud.Count);
% atlas needs to match volume dimensions
% tv    = permute(tv, [1 3 2]);

% tvaffine = imwarp(tv, Rmoving, tform, 'OutputView',Rfixed);
%==========================================================================
% the first step is to make sure our sample is nicely aligned
fprintf('Obtaining initial similarity transform... '); tic;
transinit = originalSimilarityTform(ls_cloud, tv_cloud, params, Tflip);
fprintf('Done! Took %2.1f s. \n', toc);

fprintf('Identifying candidate corresponding points... '); tic;
[cpsample, cpatlas] = triageAndMatchClouds(ls_cloud, tv_cloud, transinit, params);
fprintf('Done! Took %2.1f s. \n', toc);
%==========================================================================
% we plot the first step
Rtemplate = imref3d(size(tvreg));
Rsample   = imref3d(size(volumereg));
avtest    = imwarp(avreg, Rtemplate, transinit.invert,'nearest', 'OutputView',Rsample);
volmax    = quantile(volumereg,0.999,'all');
for idim = 1:3
    cf = plotAnnotationComparison(uint8(255*volumereg/volmax), avtest, idim);
    print(cf, fullfile(opts.savepath, sprintf('dim%d_initial_registration', idim)), '-dpng')
    close(cf);
end
%==========================================================================
% let's save stuff
opts.permute_sample_to_atlas = [1 3 2];
opts.original_trans          = transinit;
opts.downfac_reg             = downfac;
opts.autocpsample            = cpsample;
opts.autocpatlas             = cpatlas;
% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')
%==========================================================================
end