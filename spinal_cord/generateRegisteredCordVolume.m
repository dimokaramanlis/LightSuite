function straightvol = generateRegisteredCordVolume(regopts, transformparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
cordvol   = readSpinalCordSample(regopts.datafolder, regopts.sampleres);
Nchannels = size(cordvol, 4);
%==========================================================================
channames = cell(Nchannels, 1);
if isfield(regopts, 'channames')
    channames = regopts.channames;
end
%==========================================================================
% we prepare the cord volume
howtoperm = [transformparams.how_to_perm 4];
yrange    = transformparams.samp_ikeepy;
xrange    = transformparams.samp_ikeepx;
zrange    = transformparams.samp_ikeeplong;

cordvol   = permute(cordvol, howtoperm);
cordvol   = cordvol(yrange(1):yrange(2), xrange(1):xrange(2), zrange(1):zrange(2), :);
%==========================================================================
fprintf('Applying transforms... \n'); savetic = tic;

tforms        = transformparams.slicetforms;
bsplineparams = transformparams.tform_bspline_samp20um_to_atlas_20um_px;
affineparams  = transformparams.tform_affine_samp20um_to_atlas_20um_px;

sizetv     = transformparams.atlassize([1 2]);
raout      = imref2d(sizetv);
sizesample = [sizetv numel(tforms)];
Rmoving    = imref3d(sizesample);
Rfixed     = imref3d(transformparams.atlassize);

straightvol = zeros(transformparams.atlassize, 'uint16');
for ichan = 1:Nchannels
    currvol    = tranformCordImagesSlices(cordvol(:, :, :, ichan), tforms, raout);
    volumereg  = transformix(currvol, bsplineparams,...
    'movingscale', 0.02*[1 1 1]);
    straightvol(:, :, :, ichan) = ...
        imwarp(volumereg, Rmoving, affineparams, 'OutputView', Rfixed);
    fprintf('Channel %d/%d done. Took %2.2f s. \n', ichan, Nchannels, toc(savetic));
end
%==========================================================================
% we finally flip along the long axis if needed
if transformparams.tofliprc
    straightvol = flip(straightvol, 3);
end
%==========================================================================
fprintf('Saving registered volume... '); savetic = tic;
finalpath     = fullfile(regopts.lsfolder, 'volume_registered');
saveLargeSliceVolume(permute(straightvol, [1 2 4 3]), channames, finalpath);
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