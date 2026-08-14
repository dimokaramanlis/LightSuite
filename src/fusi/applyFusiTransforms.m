function [res] = applyFusiTransforms(opts, data, voxelsize_mm, savefullvols, transfuntoanatomy)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
if nargin < 5 || isempty(transfuntoanatomy)
    warning('applyFusiTransforms:noRigidTransform', ...
        'No functional->anatomy transform given; assuming data is already in anatomy space')
    rigidfirst = false;
else
    rigidfirst = true;
end
%--------------------------------------------------------------------------
Nframes         = size(data, 4);
transformparams = fullfile(opts.savepath,"transform_params.mat");
trstruct        = load(transformparams);
permvec         = opts.permute_sample_to_atlas;
facmovie        = voxelsize_mm(abs(permvec))./opts.atlas_res;
data            = permuteBrainVolume(data, permvec);
movievolsize    = ceil(facmovie.*size(data,1:3));
Rmovie          = imref3d(movievolsize);
%--------------------------------------------------------------------------
% prepare non-rigid transform
transpath = trstruct.tform_bspline_samp20um_to_atlas_20um_px;
paramsfin = elastix_parameter_read(transpath);
paramsfin.ResultImagePixelType = 'float';
[savepath, fname, fext] = fileparts(transpath);
new_temp_path = fullfile(savepath, sprintf('%s_temp_annotation%s', fname, fext));
elastix_paramStruct2txt(new_temp_path, paramsfin);
%--------------------------------------------------------------------------
ratalas    = imref3d(trstruct.atlassize);

facanatomy = opts.sample_res(abs(permvec))./opts.atlas_res;
rasamp     = size(opts.sample, 1:3);
rasamp     = rasamp(abs(permvec));
rasamp     = imref3d(ceil(rasamp.*facanatomy));
% 
groupidx   = unique(opts.parcelinfo.parcellation_index);
Nforaccum  = max(opts.annotation, [], 'all') + 1;

res.areasignals    = zeros(Nforaccum, Nframes, 'single'); % also hemisphere
res.areasignalsaff = zeros(Nforaccum, Nframes, 'single'); % also hemisphere
areavols           = accumarray(opts.annotation(:)+1, 1, [Nforaccum 1], @sum);
areavols           = areavols * 0.01^3; % in mm
% [~, ~, zall]       = meshgrid(1:trstruct.atlassize(1), 1:trstruct.atlassize(2), 1:trstruct.atlassize(3));
% sideid             = (zall>=trstruct.atlassize(3)/2);
% 

msg = [];    tic;

if savefullvols
    res.vreg    = nan([trstruct.atlassize, Nframes], 'single');
    res.vregaff = nan([trstruct.atlassize, Nframes], 'single');
end

% -------------------------------------------------------------------------
% Precompute the dense B-spline displacement field ONCE. It only depends on
% the transform, not on the data, so we pay the transformix/elastix cost a
% single time and reuse the field with imwarp for every frame.
%
% transformix -def all (empty moving image) returns the elastix deformation
% field u(x) = T(x) - x for each voxel of the fixed grid, in PHYSICAL units
% (mm). mhd_read lays a vector field out as [component, x, y, z], with the
% vector components AND the spatial axes both in elastix (x, y, z) order (no
% [2 1 3] row/col swap is applied to vector fields, unlike scalar images).
%
% imwarp wants a field of size [rows cols planes 3] = [y x z comp] in
% INTRINSIC (voxel) units, component order (1=x=col, 2=y=row, 3=z=plane),
% used as a backward / "pull" map:  sample_location = voxel_location + D.
% elastix resamples with the same pull convention, Result(x)=Moving(x+u(x)),
% so there is NO sign flip (the earlier -D was the bug). Two conversions:
%   1) mm -> voxels: divide each component by the fixed-grid spacing.
%   2) reorder [comp x y z] -> [y x z comp] via permute([3 2 4 1]); this
%      swaps the x/y spatial axes (elastix col-major vs MATLAB row-major) and
%      moves the component to dim 4, leaving comp order (x,y,z) intact.
% -------------------------------------------------------------------------
Dfield = transformix([], new_temp_path, 'verbose', false);  % [3 x y z], mm
Dfield = Dfield ./ paramsfin.Spacing(:);                    % mm -> voxels (per comp)
Dfield = permute(Dfield, [3 2 4 1]);                        % [y x z comp], voxels
%-------------------------------------------------------------------------
% for fast accumulation
outcell = accumarray([opts.annotation(:)]+1, ...
    1:prod(trstruct.atlassize), [Nforaccum 1], @(x) {x});
%-------------------------------------------------------------------------

for iframe = 1:Nframes
    %----------------------------------------------------------------------
    % prepare frame
    currvol   = data(:, :, :, iframe);
    currvol   = imresize3(currvol, 'Scale', facmovie);
    %----------------------------------------------------------------------
    if rigidfirst
        currvol   = imwarp(currvol, Rmovie, transfuntoanatomy.tformrigid,...
            'OutputView', transfuntoanatomy.Ranatomy, 'FillValues',nan);
    end
    %----------------------------------------------------------------------
    % non-rigid B-spline, applied via the precomputed displacement field.
    % 'linear' is fast and robust to the NaN FOV border; switch to 'cubic'
    % to match transformix's FinalBSplineInterpolationOrder = 3 more closely.
    volumereg = imwarp(currvol, Dfield, 'cubic', 'FillValues', nan);

    % volumereg = transformix(currvol,new_temp_path,...
    %     'movingscale', opts.atlas_res.*[1 1 1], 'verbose', false);
    %----------------------------------------------------------------------
    volumereg = imwarp(volumereg, rasamp, trstruct.tform_affine_samp20um_to_atlas_10um_px, ...
        'OutputView',ratalas, 'FillValues',nan);
    %----------------------------------------------------------------------
    % affine only
    volumeregaff = imwarp(currvol, rasamp, trstruct.tform_affine_samp20um_to_atlas_10um_px, ...
        'OutputView',ratalas, 'FillValues',nan);
    %----------------------------------------------------------------------
    % res.areasignals(:,    iframe) = accumarray([opts.annotation(:)]+1, ...
    %     volumereg(:),    [Nforaccum 1], @nanmean);
    % res.areasignalsaff(:, iframe) = accumarray([opts.annotation(:)]+1, ...
    %     volumeregaff(:), [Nforaccum 1], @nanmean);

    res.areasignals(:,    iframe) = cellfun(@(x) mean(volumereg(x), 'omitmissing'), outcell);
    res.areasignalsaff(:, iframe) = cellfun(@(x) mean(volumeregaff(x), 'omitmissing'), outcell);
    %----------------------------------------------------------------------
    if savefullvols
        res.vreg(:, :, :, iframe)    = volumereg;
        res.vregaff(:, :, :, iframe) = volumeregaff;
    end
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Frame %d/%d. Time elapsed %2.2f s...\n',...
        iframe, Nframes,toc);
    fprintf(msg);
    %----------------------------------------------------------------------
end
delete(new_temp_path);
%----------------------------------------------------------------------
res.areasignals    = res.areasignals(groupidx+1, :);
res.areasignalsaff = res.areasignalsaff(groupidx+1, :);
res.areavols       = areavols(groupidx+1);
res.groupidx       = groupidx;
%==========================================================================
end