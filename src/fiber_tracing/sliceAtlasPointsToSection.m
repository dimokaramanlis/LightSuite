function [isec, rc, Dfields] = sliceAtlasPointsToSection(pts, trstruct, sliceinfo, savepath, Dfields)
% SLICEATLASPOINTSTOSECTION  Map Allen CCF points to pixels of the aligned
%   sections (volume_aligned/, px_process) of a slice-module brain.
%
%   [ISEC, RC, DFIELDS] = sliceAtlasPointsToSection(PTS, TRSTRUCT, SLICEINFO, SAVEPATH, DFIELDS)
%
%   Follows generateRegisteredSliceVolume pixel for pixel, so a point lands
%   on the section pixel that volume_registered/ shows there:
%     atlas 10 µm voxel -> 20 µm atlas world -> rigid (allen_to_samp) ->
%     sample-space plane (nearest section, as fillmissing 'nearest') ->
%     per-section affine (atlas_to_image) -> inverse b-spline -> section px.
%
%   PTS       N×3 [AP DV ML] full-CCF 10 µm voxel coordinates (fractional ok).
%   TRSTRUCT  transform_params.mat of the slice brain.
%   SLICEINFO sliceinfo struct (px_process, px_register, px_atlas).
%   DFIELDS   optional cache of the per-section b-spline displacement fields
%             (cell, filled on demand and returned).
%
%   ISEC  N×1 section index (NaN outside the first/last section).
%   RC    N×2 [row col] in the volume_aligned/ section image (1-based).

    Nsec = numel(trstruct.sliceids_in_sample_space);
    if nargin < 5 || isempty(Dfields); Dfields = cell(1, Nsec); end

    % registered-volume grid: the 20 µm atlas reference at 10 µm pixels
    Ra  = trstruct.space3d_atlas_20um;
    Ra  = imref3d(Ra.ImageSize * 2, Ra.XWorldLimits, Ra.YWorldLimits, Ra.ZWorldLimits);
    Rs0 = trstruct.space3d_sample_20um;
    sz2 = Ra.ImageSize(2:3);                          % in-plane [DV ML] at 10 µm
    Rs  = imref3d([ceil(Rs0.ImageExtentInWorldY), sz2], ...
        Rs0.XWorldLimits, Rs0.YWorldLimits, Rs0.ZWorldLimits);

    % atlas voxel -> atlas world (imref3d: rows = Y = AP, cols = X = DV)
    apr = pts(:, 1) - (trstruct.atlasaplims(1) - 1);
    [xw, yw, zw] = intrinsicToWorld(Ra, pts(:, 2), apr, pts(:, 3));
    % imwarp(..., rigid.invert) samples the sample volume at rigid(w)
    sw = trstruct.tformrigid_allen_to_samp_20um.transformPointsForward([xw yw zw]);
    [cdv, rap, pml] = worldToIntrinsic(Rs, sw(:, 1), sw(:, 2), sw(:, 3));

    % nearest section along the sample AP axis; none beyond the end sections
    sids = trstruct.sliceids_in_sample_space(:)';
    [~, isec] = min(abs(rap - sids), [], 2);
    isec = double(isec);
    isec(rap < min(sids) - 0.5 | rap > max(sids) + 0.5) = NaN;

    rc = nan(size(pts, 1), 2);
    sc = sliceinfo.px_register / sliceinfo.px_atlas;   % 10 µm px per 20 µm unit
    for k = unique(isec(~isnan(isec)))'
        ii = isec == k;
        % plane pixel (row = DV, col = ML) -> rahist world -> affine
        w  = [pml(ii), cdv(ii)] / sc;
        u  = trstruct.tformaffine_tform_atlas_to_image(k).transformPointsForward(w);
        b  = u * sc;                                   % [col row] on the 10 µm grid
        if isempty(Dfields{k})
            Dfields{k} = sectionDisplacement(trstruct, sliceinfo, savepath, k, sz2);
        end
        D  = Dfields{k};                               % [2 x Ncol x Nrow], mm
        dx = interpn(D{1}, b(:, 1), b(:, 2), 'linear', 0);
        dy = interpn(D{2}, b(:, 1), b(:, 2), 'linear', 0);
        % elastix grid: origin 0, spacing px_atlas; moving image spacing px_process
        m  = ((b - 1) * sliceinfo.px_atlas * 1e-3 + [dx dy]) / (sliceinfo.px_process * 1e-3) + 1;
        rc(ii, :) = m(:, [2 1]);
    end
end

function D = sectionDisplacement(trstruct, sliceinfo, savepath, k, sz2)
% b-spline displacement of section K on the 10 µm atlas plane, exactly the
% grid generateRegisteredSliceVolume resamples on
    revdir = fullfile(savepath, 'elastix_reverse');
    if ~exist(revdir, 'dir'); revdir = trstruct.tformbspline_samp20um_to_atlas_20um; end
    f  = dir(fullfile(revdir, sprintf('%03d_slice_*', k)));
    el = elastix_parameter_read(fullfile(f(1).folder, f(1).name));
    el.Size    = flip(sz2);
    el.Spacing = ones(1, 2) * sliceinfo.px_atlas * 1e-3;
    tmp = [tempname '.txt'];
    evalc('elastix_paramStruct2txt(tmp, el)');       % silence its log line
    cl  = onCleanup(@() delete(tmp));
    F   = transformix([], tmp, 'verbose', false);       % [2 x Ncol x Nrow]
    D   = {squeeze(F(1, :, :)), squeeze(F(2, :, :))};
end
