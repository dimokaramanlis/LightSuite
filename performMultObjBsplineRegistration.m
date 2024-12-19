function [regimg,tform_bspline, tformpath, pathtemp] = performMultObjBsplineRegistration(movingvol,fixedvol,volscale,...
    movingpts, fixedpts, cpwt, savepath)
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
params.FixedImageDimension             = 3;
params.MovingImageDimension            = 3;
params.FixedImagePyramid               = 'FixedRecursiveImagePyramid';
params.MovingImagePyramid              = 'MovingRecursiveImagePyramid';
params.UseRandomSampleRegion           = true;
params.NewSamplesEveryIteration        = true;
params.NumberOfResolutions             = 4;
params.NumberOfHistogramBins           = 32;
params.SP_A                            = 2000;
%--------------------------------------------------------------------------
% these may affect more
params.MaximumNumberOfIterations       = [600   600 1200 1200]; % 1000 was also good
params.NumberOfSpatialSamples          = [1000 1000 2000 2000];%[500 1000 2000 4000]; % 3e3
params.Metric1Weight                   = cpwt;
params.Metric0Weight                   = 1-cpwt;
params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];
params.FinalGridSpacingInPhysicalUnits = 0.5*ones(1,3); % 45 is consistently good
params.SampleRegionSize                = 2.5*ones(1,3);%2.5*ones(1,3); 4 was good
%--------------------------------------------------------------------------

pathtemp = fullfile(savepath, 'elastix_temp');
makeNewDir(pathtemp);

movpath  = fullfile(pathtemp, 'moving.txt');
fixpath  = fullfile(pathtemp, 'fixed.txt');

% -1 because offset is always zero. this way, the first pixel is at 0
writePointsFile(movpath, (movingpts-1)*volscale)
writePointsFile(fixpath, (fixedpts-1)*volscale)


fprintf('Performing B-spline registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml','paramstruct',params,...
    'movingpoints', movpath,  'fixedpoints', fixpath,...
    'movingscale', volscale*[1 1 1],  'fixedscale', volscale*[1 1 1]);
fprintf('Done! Took %2.2f s.\n', toc)

% copy file and delete temporary folder
tformpath = fullfile(savepath, 'bspline_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
% rmdir(pathtemp, 's');
%==========================================================================
end