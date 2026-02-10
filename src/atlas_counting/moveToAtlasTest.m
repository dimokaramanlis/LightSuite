function medianoverareas = moveToAtlasTest(inputpts, inputvol, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------

allen_atlas_path_ori = fileparts(which('template_volume_10um.npy'));
% avori               = readNPY(fullfile(allen_atlas_path_ori,'annotation_volume_10um_by_index.npy'));
tv     = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
%--------------------------------------------------------------------------
% we permute the volume to match atlas
volume  = permute(inputvol, trstruct.how_to_perm);
volumereg = transformix(volume,trstruct.tform_bspline_samp20um_to_atlas_20um_px,...
    'movingscale', 0.02*[1 1 1]);
volumereg          = uint16(abs(volumereg));
Rmoving            = imref3d(size(volumereg));
Rfixed             = imref3d(trstruct.atlassize);
registeredvolume   = imwarp(volumereg, Rmoving, trstruct.tform_affine_samp20um_to_atlas_10um_px,...
    'OutputView',Rfixed);
%--------------------------------------------------------------------------


outpts   = inputpts(:, 1:3) .* trstruct.ori_pxsize *1e-3;

outpts   = outpts(:, [2 1 3]);
outpts   = outpts(:, trstruct.how_to_perm);
outpts   = outpts(:, [2 1 3]);
outpts   = outpts/0.02;

Dfield = transformix([],trstruct.tform_bspline_samp20um_to_atlas_20um_px);
Dfield = permute(Dfield,[2 3 4 1]);

%%
[Sx, Sy, Sz, ~] = size(Dfield);

% Create grid vectors representing the indices (1-based)
Xgv = 1:Sx;
Ygv = 1:Sy;
Zgv = 1:Sz;

% Extract displacement components
Dx = Dfield(:,:,:,1);
Dy = Dfield(:,:,:,2);
Dz = Dfield(:,:,:,3);

dx_interpolated = interpn(Xgv, Ygv, Zgv, Dx, outpts(:,1), outpts(:,2), outpts(:,3), 'linear');
dy_interpolated = interpn(Xgv, Ygv, Zgv, Dy, outpts(:,1), outpts(:,2), outpts(:,3), 'linear');
dz_interpolated = interpn(Xgv, Ygv, Zgv, Dz, outpts(:,1), outpts(:,2), outpts(:,3), 'linear');
nan_indices = isnan(dx_interpolated) | isnan(dy_interpolated) | isnan(dz_interpolated);
if any(nan_indices)
    warning('%d points were outside the deformation field grid. Their displacement is set to zero.', sum(nan_indices));
    dx_interpolated(nan_indices) = 0;
    dy_interpolated(nan_indices) = 0;
    dz_interpolated(nan_indices) = 0;
end
interpolated_displacements = [-dx_interpolated, -dy_interpolated, -dz_interpolated]/0.02;
outpts2 = outpts + interpolated_displacements;
%%
% outpts   = (inputpts(:, 1:3) - 0.5) .* trstruct.ori_pxsize./(trstruct.atlasres*2)+0.5;
% outpts   = (inputpts(:, 1:3) - 0.5).*trstruct.regvolsize./trstruct.ori_size + 0.5;
% outpts   = inputpts(:, 1:3) .*trstruct.regvolsize./trstruct.ori_size;

outpts   = inputpts(:, 1:3) .* trstruct.ori_pxsize *1e-3;

outpts   = outpts(:, [2 1 3]);
outpts   = outpts(:, trstruct.how_to_perm);
outpts   = outpts(:, [2 1 3]);

pointlocs =  inputpts(:, 1:3) .*trstruct.regvolsize./trstruct.ori_size;
pointlocs   = pointlocs(:, [2 1 3]);
pointlocs   = round(pointlocs(:, trstruct.how_to_perm));
irem        = any(pointlocs <1,2);
pointlocs(irem,:) = [];
avplot       = ndSparse.build(double(pointlocs), 1,trstruct.regvolsize(trstruct.how_to_perm));
avplot       = single(full(avplot));
avplot       = imgaussfilt3(avplot, 1);

multfac = 2^15-1/max(avplot,[],'all');

cellsreg = transformix(uint16(avplot*multfac),trstruct.tform_bspline_samp20um_to_atlas_20um_px,...
    'movingscale', 0.02*[1 1 1]);

transout = transformix(outpts(:,[2 1 3]), trstruct.tform_bspline_samp20um_to_atlas_20um_px);

% transout = transformix(outpts(iplot,:)*0.02, trstruct.tform_bspline_samp20um_to_atlas_20um_px);
outpts2  = transout.OutputPoint(:,[2 1 3]);
outpts2  = outpts2/0.02;
outpts   = outpts/0.02;
%--------------------------------------------------------------------------
%%
clf;
islice =450;
% iplot = (outpts2(:,2) > (islice - 1)) & (outpts2(:,2) < (islice + 1));
iplot = (outpts2(:,3) > (islice - 1)) & (outpts2(:,3) < (islice + 1));

subplot(1,2,1)
imagesc(squeeze(volumereg(:,:,islice)), [100 15000]); hold on;
plot(outpts2(iplot, 1), outpts2(iplot,2), 'r.')
% imagesc(squeeze(volumereg(islice,:,:)), [100 15000]); hold on;
% plot(outpts2(iplot, 3), outpts2(iplot,1), 'r.')
subplot(1,2,2)
imagesc(squeeze(volume(:,:,islice)), [100 15000]); hold on;
plot(outpts(iplot, 1), outpts(iplot,2), 'r.')
% imagesc(squeeze(volume(islice,:,:)), [100 15000]); hold on;
% plot(outpts(iplot, 3), outpts(iplot,1), 'r.')


%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
end