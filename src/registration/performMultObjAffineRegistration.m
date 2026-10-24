function [regimg,tform_bspline, tformpath, pathtemp]  = performMultObjAffineRegistration(...
    movingvol,fixedvol,volscale,...
    movingpts, fixedpts, cpwt, savepath, paramoverrides)
%PERFORMMULTOBJAFFINEREGISTRATION affine elastix fit using image + points.
%   Works for 2D images and 3D volumes (dimension taken from FIXEDVOL).
%   Mutual information on the images is combined with the control point
%   distance, weighted by CPWT. Without control points, only the image is
%   used. PARAMOVERRIDES (struct) replaces any of the default parameters.
%   The transform is saved as SAVEPATH/affine_atlas_to_samp_20um.txt.
if nargin < 8 || isempty(paramoverrides), paramoverrides = struct(); end
%==========================================================================
addElastixRepoPaths;
params = struct();
ndim   = ndims(fixedvol);
haspts = ~isempty(movingpts) && ~isempty(fixedpts);
%==========================================================================
if haspts
    params.Registration              = 'MultiMetricMultiResolutionRegistration';
    params.Metric                    = {'AdvancedMattesMutualInformation',...
        'CorrespondingPointsEuclideanDistanceMetric'};
    params.Metric1Weight             = cpwt;
    params.Metric0Weight             = 1;
else
    params.Registration              = 'MultiResolutionRegistration';
    params.Metric                    = 'AdvancedMattesMutualInformation';
end

params.Transform                       = 'AffineTransform';%'RecursiveBSplineTransform';
params.Optimizer                       = 'AdaptiveStochasticGradientDescent';%,AdaptiveStochasticGradientDescent, QuasiNewtonLBFGSOptimizer]
params.ImageSampler                    = 'RandomCoordinate'; % RandomCoordinate
params.AutomaticParameterEstimation    = true;
params.AutomaticScalesEstimation       = true;
params.BSplineInterpolationOrder       = 3;
params.FinalBSplineInterpolationOrder  = 3;
params.FixedImageDimension             = ndim;
params.MovingImageDimension            = ndim;
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
params.ImagePyramidSchedule            = kron([8 4 2 1], ones(1, ndim));
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

elastixargs = {'movingscale', volscale*ones(1, ndim), 'fixedscale', volscale*ones(1, ndim)};
if haspts
    movpath = fullfile(pathtemp, 'moving.txt');
    fixpath = fullfile(pathtemp, 'fixed.txt');
    % -1 because offset is always zero. this way, the first pixel is at 0
    writePointsFile(movpath, (movingpts-1)*volscale)
    writePointsFile(fixpath, (fixedpts-1)*volscale)
    elastixargs = [elastixargs, {'movingpoints', movpath, 'fixedpoints', fixpath}];
end

fprintf('Performing affine registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml',...
    'paramstruct',params, elastixargs{:});
fprintf('Done! Took %2.2f s.\n', toc)
%-------------------------------------------------------------------------
% copy file and delete temporary folder
tformpath = fullfile(savepath, 'affine_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
rmdir(pathtemp, 's');
%-------------------------------------------------------------------------
end
