function opts = initializeRegistration(inputpath, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function re-orients the sample to match the atlas and facilitate
% control point selection
%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) isstring(x) || ischar(x));
addParameter(p, 'Volume', [], @isnumeric);
addParameter(p, 'outlierratio', 0.01, @isnumeric);
addParameter(p, 'cloudthres', 5, @isnumeric);

parse(p, inputpath, varargin{:});
params = p.Results;
%==========================================================================
dpbcpd          = which('bcpd.exe');
params.bpcdpath = dpbcpd;
%==========================================================================
optsfile   = dir(fullfile(inputpath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
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
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
tvreg            = imresize3(tv, downfac);
%==========================================================================
bofile = fullfile(opts.savepath, 'brain_orientation.txt');
fprintf('Looking for brain orientation data in %s\n', bofile)
if exist(bofile, 'file')
    permvec = load(bofile);
    fprintf('Found it, brain orientation is %s\n', mat2str(permvec))
else
    fprintf('You have to specify the orientation, check GUI\n')
    permvec = getBrainOrientation(backvol,tvreg);
    % Save to file
    writematrix(permvec, bofile);
    fprintf('Brain orientation saved as %s to %s\n', mat2str(permvec), bofile);
end
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
fprintf('Creating cloud for sample volume... '); tic;
volumereg  = permuteBrainVolume(newvol, permvec);
ls_cloud   = extractSamplePoints(volumereg, params.cloudthres);
fprintf('Done! Took %2.1f s. Found %d points.\n', toc, ls_cloud.Count);
%==========================================================================
% load atlas and extract corresponding points
fprintf('Loading atlas data and generating the atlas cloud... '); tic;
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
avreg   = imresize3(av, downfac, 'Method','nearest');

% tv_cloud = extractVolumePoints(tvreg, 15);
tvforpoints             = single(tvreg);
tvforpoints(avreg == 0) = 0;
tv_cloud                = extractVolumePointsGradient(tvforpoints, 20, 5);
fprintf('Done! Took %2.1f s. Found %d points.\n', toc, tv_cloud.Count);
%==========================================================================
% the first step is to make sure our sample is nicely aligned
fprintf('Obtaining initial similarity transform... '); tic;
transinit = originalSimilarityTform(ls_cloud, tv_cloud, params);
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
opts.permute_sample_to_atlas = permvec;
opts.original_trans          = transinit;
opts.downfac_reg             = downfac;
opts.autocpsample            = cpsample;
opts.autocpatlas             = cpatlas;
% save registration
save(fullfile(opts.savepath, 'regopts.mat'), 'opts')
%==========================================================================
end


