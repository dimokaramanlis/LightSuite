function [regimg,tform_bspline] = bsplineRegisterSlice(movingim,fixedim,...
    imscale, movingpts, fixedpts, cpwt, savepath, gridspacing, verbose)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   BSPLINEREGISTERSLICE(..., GRIDSPACING) sets the final B-spline control
%   point spacing in mm (how much the slice is allowed to bend); smaller
%   values allow more local deformation. Defaults to 0.96 when not given.
%   A scalar is used for both axes; [sx sy] sets x and y separately.
%   VERBOSE (default true) prints progress and timing.
%==========================================================================
if nargin < 8 || isempty(gridspacing)
    gridspacing = 0.96;
end
if nargin < 9 || isempty(verbose), verbose = true; end
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
params.NumberOfHistogramBins           = 32;
params.SP_A                            = 20;
%--------------------------------------------------------------------------
% these may affect more
params.FinalGridSpacingInPhysicalUnits = gridspacing(:)' .* ones(1,2); % scalar or [x y]
params.MaximumNumberOfIterations       = [600 800 1200 1600];
params.NumberOfSpatialSamples          = 5000;
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


if verbose, fprintf('Performing B-spline registration with elastix...\n'); end
tic;
[regimg,tform_bspline] = elastix(movingim, fixedim, savepath,'elastix_default.yml','paramstruct',params,...
    'movingpoints', movpath,  'fixedpoints', fixpath,...
    'movingscale', imscale*[1 1],  'fixedscale', imscale*[1 1]);
if verbose, fprintf('Done! Took %2.2f s.\n', toc); end
%==========================================================================
end