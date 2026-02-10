function [regimg,tform_bspline, tformpath, pathtemp]  = performMultObjAffineRegistration(movingvol,fixedvol,volscale,...
    movingpts, fixedpts, cpwt, savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
addElastixRepoPaths;
params = struct();
%==========================================================================
params.Registration                  = 'MultiMetricMultiResolutionRegistration';
params.Metric                        = {'AdvancedMattesMutualInformation',...
    'CorrespondingPointsEuclideanDistanceMetric'};
params.Metric1Weight                   = cpwt;
params.Metric0Weight                   = 1;

params.Transform                       = 'AffineTransform';%'RecursiveBSplineTransform';
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
params.NumberOfHistogramBins           = 16;
%--------------------------------------------------------------------------
% these may affect more
params.MaximumNumberOfIterations       = [500 1000 1500 2000]; %1000; %[1000 1500 2000 2500]; %
params.NumberOfSpatialSamples          = 5000;%[1000 1000 2000 2000];% [2000 2500 3000 3000];%
params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];
%--------------------------------------------------------------------------
pathtemp = fullfile(savepath, 'elastix_temp');
makeNewDir(pathtemp);
movpath = fullfile(pathtemp, 'moving.txt');
fixpath = fullfile(pathtemp, 'fixed.txt');


writePointsFile(movpath, (movingpts-1)*volscale)
writePointsFile(fixpath, (fixedpts-1)*volscale)


fprintf('Performing affine registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml','paramstruct',params,...
    'movingpoints', movpath,  'fixedpoints', fixpath,...
    'movingscale', volscale*[1 1 1],  'fixedscale', volscale*[1 1 1]);
fprintf('Done! Took %2.2f s.\n', toc)
%-------------------------------------------------------------------------
% copy file and delete temporary folder
tformpath = fullfile(savepath, 'affine_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
rmdir(pathtemp, 's');
%-------------------------------------------------------------------------
end