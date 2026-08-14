function [regimg,tform_bspline, tformpath, pathtemp] = performElastixRigidRegistration(movingvol,fixedvol,...
    volscale, savepath, paramoverrides)
%PERFORMELASTIXRIGIDREGISTRATION Rigid (Euler) registration with elastix.
%   paramoverrides (optional struct) overwrites individual elastix parameters
%   set below, e.g. struct('NumberOfResolutions',1,'ImagePyramidSchedule',[1 1 1])
%   to fit at a single full scale when the inputs are already roughly aligned.
%==========================================================================
if nargin < 5 || isempty(paramoverrides), paramoverrides = struct(); end
addElastixRepoPaths;
params = struct();
%==========================================================================
% general parameters, probably won't touch
params.Registration                    = 'MultiResolutionRegistration';
params.Metric                          = 'AdvancedMattesMutualInformation';
params.Transform                       = 'EulerTransform';%'EulerTransform';%'RecursiveBSplineTransform';
params.Optimizer                       = 'AdaptiveStochasticGradientDescent';%,AdaptiveStochasticGradientDescent, QuasiNewtonLBFGSOptimizer]
params.ImageSampler                    = 'RandomCoordinate'; % RandomCoordinate
params.AutomaticParameterEstimation    = true;
params.AutomaticScalesEstimation       = true;
params.BSplineInterpolationOrder       = 3;
params.FinalBSplineInterpolationOrder  = 3;
params.FixedImageDimension             = 3;
params.MovingImageDimension            = 3;
params.FixedImagePyramid               = 'FixedRecursiveImagePyramid';
params.MovingImagePyramid              = 'MovingRecursiveImagePyramid';
params.UseRandomSampleRegion           = true;
params.NewSamplesEveryIteration        = true;
params.NumberOfResolutions             = 4;
params.NumberOfHistogramBins           = 32;
%--------------------------------------------------------------------------
% these may affect more
params.MaximumNumberOfIterations       = [500 1000 1500 2000]; %1000; %[1000 1500 2000 2500]; %
params.NumberOfSpatialSamples          = 5000;%[1000 1000 2000 2000];% [2000 2500 3000 3000];%
params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];
%--------------------------------------------------------------------------
% caller overrides (e.g. a single full-scale resolution for a warm-started
% residual fit) take precedence over the defaults above
fn = fieldnames(paramoverrides);
for ii = 1:numel(fn)
    params.(fn{ii}) = paramoverrides.(fn{ii});
end
%--------------------------------------------------------------------------

pathtemp = fullfile(savepath, 'elastix_temp');
makeNewDir(pathtemp);

fprintf('Performing rigid registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml','paramstruct',params,...
    'movingscale', volscale.*[1 1 1],  'fixedscale', volscale.*[1 1 1]);
fprintf('Done! Took %2.2f s.\n', toc)

% copy file and delete temporary folder
tformpath = fullfile(savepath, 'rigid_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
rmdir(pathtemp, 's');
%==========================================================================
end