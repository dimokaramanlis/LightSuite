function tformsout = alignConsecutiveCordSlicesImage(regopts, tformslices)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


Nslices = numel(tformslices);
%==========================================================================
sizetv   = size(regopts.tv, [2 3]);
tvxlims  = sizetv(2)/2 + [-1 1]*sizetv(2);
tvylims  = sizetv(1)/2 + [-1 1]*sizetv(1);
raout    = imref2d([range(tvylims) range(tvxlims)], tvxlims, tvylims);
sampletr = tranformCordImagesSlices(regopts.regvol, tformslices, raout);
%==========================================================================
Nsamples   = min(100, Nslices);
minstart   = Nsamples/2;
maxend     = Nslices - Nsamples/2;
idxsamples = unique(round(linspace(minstart+1, maxend-1, Nsamples)));

tformarray(numel(idxsamples),1) = rigidtform2d; % initialize as identity
tformarray(1) = rigidtform2d;
for isamp = 1:numel(idxsamples)-1
    slicetic = tic;
    fprintf('Aligning batch %d (moving) to %d (fixed)... ', isamp+1, isamp); 

    fixedim  = squeeze(sampletr(idxsamples(isamp), :, :));
    movingim = squeeze(sampletr(idxsamples(isamp+1), :, :));

    [optimizer,metric] = imregconfig("multimodal");
    % optimizer.MaximumIterations = 600;
    % metric.NumberOfHistogramBins = 10;
    tformout   = imregtform(movingim, fixedim, 'rigid', optimizer,metric);
    
    % newim   = imwarp(movingim, raout, tformout, 'OutputView',raout);
    % fimshow = uint8(255 * single(fixedim)/single(quantile(fixedim,0.999,'all')));
    % imshowpair(fimshow, newim)
    %  pause;
    tformarray(isamp+1) = rigidtform2d(tformarray(isamp).A * tformout.A);
    fprintf('Done! Took %2.2f s\n', toc(slicetic));

    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
end
%==========================================================================
% interpolate between batches to get back one per slice
reltrans    = reshape([tformarray.Translation], 2, numel(idxsamples))';
reltrans    = reltrans - median(reltrans);
Nfilter     = max(floor(numel(idxsamples)/6)*2+1, 3);
reltrans    = smoothdata(reltrans, 'rlowess', Nfilter);
relangles   = deg2rad([tformarray.RotationAngle])';
relvectors  = [sin(relangles) cos(relangles)];
relvectors  = smoothdata(relvectors, 'rlowess', Nfilter);
relvectors  = interp1(idxsamples, relvectors, 1:Nslices, "pchip","extrap");
reltrans    = interp1(idxsamples, reltrans, 1:Nslices, "pchip","extrap");

relangles2  = rad2deg(atan2(relvectors(:, 1), relvectors(:,2)));
tformsout   = tformslices;
for islice = 1:Nslices
    tcurr             = rigidtform2d(relangles2(islice), reltrans(islice, :));
    tformsout(islice) = rigidtform2d(tcurr.A * tformslices(islice).A);
end
%==========================================================================
end