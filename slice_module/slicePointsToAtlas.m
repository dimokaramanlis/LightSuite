function finalpts = slicePointsToAtlas(inputpts, trstruct)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
Nslices        = numel(trstruct.sliceids_in_sample_space);
bsplinedp      = trstruct.tformbspline_samp20um_to_atlas_20um;
bspltransforms = dir(fullfile(bsplinedp, '*_slice*.txt'));
bspltransforms = fullfile(bsplinedp, {bspltransforms(:).name}');
regsize_mm     = trstruct.atlasres * 2 * 1e-3;
outpts         = inputpts;
Rsample3d      = trstruct.space3d_sample_20um;
ysampvals      = linspace(Rsample3d.YWorldLimits(1)+0.5, Rsample3d.YWorldLimits(2)-0.5,...
    Rsample3d.ImageExtentInWorldY);

for islice = 1:Nslices
    %--------------------------------------------------------------------------
    currids  = inputpts(:, 3) == islice;
    currlocs = inputpts(currids, 1:3);
    pts      = currlocs(:, 1:2) * 1.25 * 1e-3/regsize_mm;
    afftrans = trstruct.tformaffine_tform_atlas_to_image(islice);
    %--------------------------------------------------------------------------
    % we first invert bsplines
    Dfield = transformix([], bspltransforms{islice});
    % we then permute and scale by the registration resolution to go to voxels
    Dfield = permute(Dfield,[2 3 1])/regsize_mm; 
    [Sx, Sy, ~] = size(Dfield);
    
    % Create grid vectors representing the indices (1-based)
    Xgv = 1:Sx;
    Ygv = 1:Sy;
    
    % Extract displacement components
    dx_interpolated = interpn(Xgv, Ygv, Dfield(:,:,1), pts(:,1), pts(:,2), 'linear');
    dy_interpolated = interpn(Xgv, Ygv, Dfield(:,:,2), pts(:,1), pts(:,2), 'linear');
    
    % negative sign is extremely important!!!
    interpolated_displacements = -[dx_interpolated, dy_interpolated];
    interpolated_displacements(isnan(interpolated_displacements)) = 0;
    finalpts = pts + interpolated_displacements;
    %--------------------------------------------------------------------------
    % we then invert affines
    finalpts = afftrans.transformPointsInverse(finalpts);
    outpts(currids, 1:2) = finalpts;
    outpts(currids,   3) = ysampvals(trstruct.sliceids_in_sample_space(islice));
    %--------------------------------------------------------------------------
end
%-------------------------------------------------------------------------
rigidfull          = trstruct.tformrigid_allen_to_samp_20um.invert;
atlasptcoords      = rigidfull.transformPointsForward(outpts(:, [2 3 1]))*2;
atlasptcoords(:,2) = atlasptcoords(:,2) + trstruct.atlasaplims(1);
%-------------------------------------------------------------------------
finalpts = [atlasptcoords inputpts(:, 4:end)];
%-------------------------------------------------------------------------
end