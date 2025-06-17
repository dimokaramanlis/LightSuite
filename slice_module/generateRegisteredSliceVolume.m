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
rahist        = imref2d(atlassizepx, 0.5, 0.5);

finvol = nan([Ny, atlassizepx, Nchannels], 'single');

for islice = 1:Nslices
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
        testim  = imwarp(testim, ratlas, slicetforms(islice).invert, 'OutputView',ratlas);
        finvol(atlasinds(islice), :, :, ichan) = testim;
    end
    %----------------------------------------------------------------------
    delete(temppath);
    %----------------------------------------------------------------------
end
%==========================================================================
% we here smooth the final volume in between slices and transform it to the
% atlas

finvol2 = smoothdata(finvol,1, 'gaussian', 7.5);
finvol2 = fillmissing(finvol, 'nearest', 1, 'EndValues','none');

regvol  = imwarp(finvol2,  Rout, opts.tformrigid_allen_to_samp_20um.invert, 'linear', 'OutputView', Ratlas);
%==========================================================================
% we finally save the registered volume

%==========================================================================
subplot(1,3,1)
imagesc(squeeze(finvol(:,120,:)))
axis equal; axis tight;
subplot(1,3,2)
imagesc(squeeze(volout(:,120,:)))
axis equal; axis tight;
subplot(1,3,3)
imagesc(squeeze(tvdown(:,120,:)))
axis equal; axis tight;

subplot(1,3,1)
imagesc(squeeze(finvol(:,:,220)),[100 2e4])
axis equal; axis tight;
subplot(1,3,2)
imagesc(squeeze(volout(:,:,220)),[100 2e4])
axis equal; axis tight;
subplot(1,3,3)
imagesc(squeeze(tvdown(:,:,220)))
axis equal; axis tight;
% % imshowpair(testim, atlasim)
%==========================================================================

end