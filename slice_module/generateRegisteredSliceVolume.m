function regvol = generateRegisteredSliceVolume(sliceinfo, transformparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
fprintf('Loading data in memory... '); tic;
slicevol             = loadLargeSliceVolume(sliceinfo.slicevolfin);
[Nchannels, Nslices] = size(slicevol, [3 4]);
fprintf('Done! Took %2.2f s\n', toc); 
%==========================================================================
% we recreate a 3d volume
atlasspacefit = transformparams.space_atlas_20um;
revsavepath   = transformparams.tformbspline_samp20um_to_atlas_20um;
slicetforms   = transformparams.tformaffine_tform_atlas_to_image;
movscale      = ones(1, 2) * sliceinfo.px_process * 1e-3;
atlassizepx   = atlasspacefit.ImageSize * sliceinfo.px_register/sliceinfo.px_atlas;
atlasinds     = transformparams.sliceids_in_sample_space;
rahist        = imref2d(atlassizepx, 0.5, 0.5);
Nysamp        = transformparams.space3d_sample_20um.ImageExtentInWorldY;

finvol = nan([Nysamp, atlassizepx, Nchannels], 'single');

for islice = 1:Nslices
    fprintf('Transforming slice %d/%d...\n', islice, Nslices);
    %----------------------------------------------------------------------
    histim           = squeeze(slicevol(:, :, :, islice));
    %----------------------------------------------------------------------
    % prepare inverse bspline
    revtformpath     = dir(fullfile(revsavepath, sprintf('%03d_slice_*', islice)));
    revtformpath     = fullfile(revtformpath.folder, revtformpath.name);
    elparams         = elastix_parameter_read(revtformpath);
    elparams.Size    = flip(atlasspacefit.ImageSize * sliceinfo.px_register/sliceinfo.px_atlas);
    elparams.Spacing = ones(1, 2) * sliceinfo.px_atlas * 1e-3;

    temppath = fullfile(sliceinfo.procpath, 'temp.txt');
    elastix_paramStruct2txt(temppath, elparams);
    %----------------------------------------------------------------------
    % apply inverse bspline and inverse affine
    for ichan = 1:Nchannels
        testim  = transformix(histim(:,:,ichan), temppath, 'movingscale',movscale, 'verbose', 0);
        testim  = uint16(abs(testim));
        testim  = imwarp(testim, rahist, slicetforms(islice).invert, 'OutputView',rahist);
        finvol(atlasinds(islice), :, :, ichan) = testim;
    end
    %----------------------------------------------------------------------
    delete(temppath);
    %----------------------------------------------------------------------
end
%==========================================================================
% smooth in between slices, transform and save

fprintf('Rigidly placing sample back to Allen... '); 
savetic = tic; msg     = [];
Routori = transformparams.space3d_sample_20um;
Ratlas  = transformparams.space3d_atlas_20um;
Ratlas  = imref3d(Ratlas.ImageSize*2, Ratlas.XWorldLimits, Ratlas.YWorldLimits, Ratlas.ZWorldLimits);
tform3d = transformparams.tformrigid_allen_to_samp_20um.invert;
regvol        = zeros([Ratlas.ImageSize([1 2]) Nchannels Ratlas.ImageSize(3)], 'uint16');
for ichan = 1:Nchannels
    % fill in gaps between slices
    finvol2 = fillmissing(finvol(:, :, :, ichan), 'nearest', 1, 'EndValues','none');
    finvol2 = uint16(finvol2);
    Rout    = imref3d(size(finvol2), Routori.XWorldLimits, Routori.YWorldLimits, Routori.ZWorldLimits);
    regvol(:, :, ichan, :)  = imwarp(finvol2,  Rout, tform3d, 'linear', 'OutputView', Ratlas);
    %------------------------------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Channel %d/%d done. Took %2.2f s. \n', ichan, Nchannels, toc(savetic));
    fprintf(msg);
    %------------------------------------------------------------------------------------------
end
%==========================================================================
fprintf('Saving registered volume... '); savetic = tic;
finalpath     = fullfile(sliceinfo.procpath, 'volume_registered');
saveLargeSliceVolume(regvol, sliceinfo.channames(ichan), finalpath);
fprintf('Done! Took %2.2f s. \n', toc(savetic));
%==========================================================================

% subplot(1,3,1)
% imagesc(squeeze(finvol(:,120,:)))
% axis equal; axis tight;
% subplot(1,3,2)
% imagesc(squeeze(volout(:,120,:)))
% axis equal; axis tight;
% subplot(1,3,3)
% imagesc(squeeze(tvdown(:,120,:)))
% axis equal; axis tight;
% 
% subplot(1,3,1)
% imagesc(squeeze(finvol(:,:,220)),[100 2e4])
% axis equal; axis tight;
% subplot(1,3,2)
% imagesc(squeeze(volout(:,:,220)),[100 2e4])
% axis equal; axis tight;
% subplot(1,3,3)
% imagesc(squeeze(tvdown(:,:,220)))
% axis equal; axis tight;
% % % imshowpair(testim, atlasim)
% %==========================================================================

end