function finalpts = volumePointsToAtlas(inputpts, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
regsize_mm = trstruct.atlasres * 2 * 1e-3;
%--------------------------------------------------------------------------
% we first transform points to match the downsampled volume
pts   = (inputpts(:, 1:3) - 1) .* trstruct.ori_pxsize *1e-3; % is this correct???

pts   = pts(:, [2 1 3]);
pts   = pts(:, trstruct.how_to_perm);
pts   = pts(:, [2 1 3]);
pts   = pts/regsize_mm;
%--------------------------------------------------------------------------
% we calculate the displacement field, which comes in world coordinates
Dfield = transformix([],trstruct.tform_bspline_samp20um_to_atlas_20um_px);
% we then permute and scale by the registration resolution to go to voxels
Dfield = permute(Dfield,[2 3 4 1])/regsize_mm; 
%--------------------------------------------------------------------------
[Sx, Sy, Sz, ~] = size(Dfield);

% Create grid vectors representing the indices (1-based)
Xgv = 1:Sx;
Ygv = 1:Sy;
Zgv = 1:Sz;

% Extract displacement components
dx_interpolated = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,1), pts(:,1), pts(:,2), pts(:,3), 'linear');
dy_interpolated = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,2), pts(:,1), pts(:,2), pts(:,3), 'linear');
dz_interpolated = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,3), pts(:,1), pts(:,2), pts(:,3), 'linear');

% negative sign is extremely important!!!
interpolated_displacements = -[dx_interpolated, dy_interpolated, dz_interpolated];
interpolated_displacements(isnan(interpolated_displacements)) = 0;
finalpts = pts + interpolated_displacements;
%--------------------------------------------------------------------------
% do we need to add anything to finalpts??? like 1 or 0.5?
finalpts = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward(finalpts);
finalpts = [finalpts inputpts(:, 4:end)];
%--------------------------------------------------------------------------
% and super finally, points should be resampled to match the 10um
% resolution of the atlas
% finalpts = 2*finalpts;
%--------------------------------------------------------------------------
% we remove points outside of the annotation volume
% irem0 = any(finalpts<0, 2);
% iremx = finalpts(:,1) > trstruct.atlassize(2);
% iremy = finalpts(:,2) > trstruct.atlassize(1);
% iremz = finalpts(:,3) > trstruct.atlassize(3);
% irem  = irem0 | iremx | iremy | iremz;
% finalpts(irem, :) = [];
%--------------------------------------------------------------------------
% nrand = min(size(finalpts,1), 1e5);
% iplot = randperm(size(finalpts,1),nrand);
% figure;
% plotBrainGrid; hold on;
% scatter3(finalpts(iplot,2),finalpts(iplot,3),finalpts(iplot,1),2,'filled','MarkerFaceAlpha',0.5)

%--------------------------------------------------------------------------
end