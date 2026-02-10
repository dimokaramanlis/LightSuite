function tformsout = alignConsecutiveCordSlices(pcsamp,  tforms, iax, batchsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
bcpd    = 'C:\Users\karamanl\Documents\GitHub\bcpd\win\bcpd.exe';
axkeep  = setdiff(1:3, iax);
Nslices        = numel(tforms);
%==========================================================================
Nbatches   = ceil(Nslices/batchsize);

tformarray(Nbatches,1) = rigidtform2d; % initialize as identity
batchcenters           = nan(Nbatches, 1);
batchcenters(1) = batchsize/2;
samppts2 = tranformCordPointsSlices(pcsamp, iax, tforms);
tformarray(1) = rigidtform2d;
for ibatch = 1:Nbatches-1
    slicetic = tic;
    fprintf('Aligning batch %d (moving) to %d (fixed)... ', ibatch+1, ibatch); 

    istart1   = (ibatch - 1) * batchsize + 1;
    iend1     = min(ibatch * batchsize, Nslices);
    ismcurr1  = samppts2(:, iax) > istart1 & samppts2(:, iax) < iend1;
    istart2   = (ibatch) * batchsize + 1;
    iend2     = min((ibatch + 1) * batchsize, Nslices);
    ismcurr2  = samppts2(:, iax) > istart2 & samppts2(:, iax) < iend2;

    pcfix   = samppts2(ismcurr1, axkeep);
    pcmov   = samppts2(ismcurr2, axkeep);
    [yreg,bfit] = pcregisterBCPD(pcmov, pcfix, 'TransformType','Rigid',...
        'BCPDPath', bcpd, 'OutlierRatio', 0.1, ...
        'Gamma', 1, Verbose = false, ConvergenceTolerance=1e-8);
    tformarray(ibatch+1) = rigidtform2d(tformarray(ibatch).A * bfit.A);
    
    batchcenters(ibatch + 1) = mean([istart2 iend2]);
    fprintf('Done! Took %2.2f s\n', toc(slicetic));
    % clf; 
    % subplot(1,2,1); hold on;
    % scatter(pcmov(:,1),pcmov(:,2),5,'k')
    % scatter(pcfix(:,1),pcfix(:,2),5, 'r')
    % axis equal; 
    % 
    % subplot(1,2,2); hold on;
    % scatter(yreg(:,1),yreg(:,2),5,'k')
    % scatter(pcfix(:,1),pcfix(:,2),5, 'r')
    % axis equal;
    % pause;
    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
end
%==========================================================================
% interpolate between batches to get back one per slice
reltrans    = reshape([tformarray.Translation], 2, Nbatches)';
Nfilter     = max(floor(Nbatches/6)*2+1, 3);
reltrans    = smoothdata(reltrans, 'rlowess', Nfilter);
relangles   = deg2rad([tformarray.RotationAngle])';
relvectors  = smoothdata([sin(relangles) cos(relangles)], 'rlowess', Nfilter);
relvectors  = interp1(batchcenters, relvectors, 1:Nslices, "pchip","extrap");
reltrans    = interp1(batchcenters, reltrans, 1:Nslices, "pchip","extrap");

relangles2  = rad2deg(atan2(relvectors(:, 1), relvectors(:,2)));
tformsout   = tforms;
for islice = 1:Nslices
    tcurr             = rigidtform2d(relangles2(islice), reltrans(islice, :));
    tformsout(islice) = rigidtform2d(tcurr.A * tforms(islice).A);
end
%==========================================================================
end