function [regimage, tform] = performMultObjAffineRegistration(movingvol,fixedvol,volscale,...
    movingpts, fixedpts, cpwt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
addElastixRepoPaths;
%==========================================================================
pathout = 'D:\lightsheet\elout';

params = struct();
params.Registration                  = 'MultiMetricMultiResolutionRegistration';
params.Metric                        = {'AdvancedMattesMutualInformation',...
    'CorrespondingPointsEuclideanDistanceMetric'};
params.Metric1Weight                   = cpwt;
params.Metric0Weight                   = 1-cpwt;
params.Transform                       = 'AffineTransform';%'RecursiveBSplineTransform';
params.FixedImageDimension             = 3;
params.MovingImageDimension            = 3;
params.NumberOfHistogramBins           = 32;
params.AutomaticTransformInitialization = true;
params.AutomaticTransformInitializationMethod = 'CenterOfGravity';
params.AutomaticScalesEstimation      = false;
% params.BSplineInterpolationOrder       = 3;
% params.FinalBSplineInterpolationOrder  = 3;
% params.Optimizer                       = 'AdaptiveStochasticGradientDescent';
% params.ImageSampler                    = 'RandomCoordinate';
% params.UseRandomSampleRegion           = true;
% params.NewSamplesEveryIteration        = true;
% params.MaximumNumberOfIterations       = [500 500 1000 1000]; % 1000 was also good
% params.NumberOfSpatialSamples          = [500 500 2000 4000];%[500 1000 2000 4000]; % 3e3

% params.NumberOfResolutions             = 4;
% params.SP_A                            = 20;
% params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];

movpath = fullfile('D:\lightsheet\elout', 'moving.txt');
fixpath = fullfile('D:\lightsheet\elout', 'fixed.txt');


writePointsFile(movpath, movingpts*volscale)
writePointsFile(fixpath, fixedpts*volscale)


tic;
[reg,stats] = elastix(movingvol, fixedvol, [],'elastix_default.yml','paramstruct',params,...
    'movingpoints', movpath,  'fixedpoints', fixpath,...
    'movingscale', volscale*[1 1 1],  'fixedscale', volscale*[1 1 1]);
toc;

delete(movpath);
delete(fixpath);
end