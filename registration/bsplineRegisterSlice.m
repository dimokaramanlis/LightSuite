function [regimg,tform_bspline] = bsplineRegisterSlice(movingim,fixedim,...
    imscale, movingpts, fixedpts, cpwt, savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
addElastixRepoPaths;
params = struct();
%==========================================================================
% general parameters, probably won't touch
params.Registration                  = 'MultiMetricMultiResolutionRegistration';
params.Metric                        = {'AdvancedMattesMutualInformation',...
    'CorrespondingPointsEuclideanDistanceMetric'};
params.Transform                       = 'RecursiveBSplineTransform';%'RecursiveBSplineTransform';
params.Optimizer                       = 'AdaptiveStochasticGradientDescent';
params.ImageSampler                    = 'RandomCoordinate';
params.AutomaticParameterEstimation    = true;
params.AutomaticScalesEstimation       = false;
params.BSplineInterpolationOrder       = 3;
params.FinalBSplineInterpolationOrder  = 3;
params.FixedImageDimension             = 2;
params.MovingImageDimension            = 2;
params.FixedImagePyramid               = 'FixedRecursiveImagePyramid';
params.MovingImagePyramid              = 'MovingRecursiveImagePyramid';
params.UseRandomSampleRegion           = true;
params.NewSamplesEveryIteration        = true;
params.NumberOfResolutions             = 4;
params.NumberOfHistogramBins           = 16;
params.SP_A                            = 20;
%--------------------------------------------------------------------------
% these may affect more
params.FinalGridSpacingInPhysicalUnits = 0.96*ones(1,2);
params.MaximumNumberOfIterations       = [600 800 1200 1600]; %1000; %[1000 1500 2000 2500]; %
params.NumberOfSpatialSamples          = 3000;%[1000 1000 2000 2000];% [2000 2500 3000 3000];%
params.Metric1Weight                   = cpwt; %cpwt;%[1  0.5 0.25 0.125] * cpwt; %
params.Metric0Weight                   = 1.0;
params.ImagePyramidSchedule            = kron([8 4 2 1],     ones(1,2));
params.SampleRegionSize                = kron([4.5 4 3 2.5], ones(1,2));
%--------------------------------------------------------------------------
makeNewDir(savepath);

movpath  = fullfile(savepath, 'moving.txt');
fixpath  = fullfile(savepath, 'fixed.txt');

% -1 because offset is always zero. this way, the first pixel is at 0
writePointsFile(movpath, (movingpts-1)*imscale)
writePointsFile(fixpath, (fixedpts-1)*imscale)


fprintf('Performing B-spline registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingim, fixedim, savepath,'elastix_default.yml','paramstruct',params,...
    'movingpoints', movpath,  'fixedpoints', fixpath,...
    'movingscale', imscale*[1 1],  'fixedscale', imscale*[1 1]);
fprintf('Done! Took %2.2f s.\n', toc)
%==========================================================================
end