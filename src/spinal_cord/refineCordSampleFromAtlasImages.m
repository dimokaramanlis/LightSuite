function [tformslices, resall] = refineCordSampleFromAtlasImages(regopts, tformslices)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Nslices = numel(tformslices);
%==========================================================================
sizetv   = size(regopts.tv, [2 3]);
tvxlims  = sizetv(2)/2 + [-1 1]*sizetv(2);
tvylims  = sizetv(1)/2 + [-1 1]*sizetv(1);
raout    = imref2d([range(tvylims) range(tvxlims)], tvxlims, tvylims);
sampletr = tranformCordImagesSlices(regopts.regvol, tformslices, raout);
%==========================================================================
atlastr  = warpCordAtlasImage(regopts.tv, regopts.atlaswarppath);
rmoving  = imref2d(sizetv);
%==========================================================================
Nsamples   = min(50, Nslices);
minstart   = Nsamples/2;
maxend     = Nslices - Nsamples/2;
idxsamples = unique(round(linspace(minstart+1, maxend-1, Nsamples)));

for isample = 1:numel(idxsamples)
    %%
    isample = 20;
    fixedim  = squeeze(sampletr(idxsamples(isample), :, :));
    movingim =  squeeze(atlastr(idxsamples(isample), :, :));
    scalefilter = 100/20;
    % [pcmov, imgmov]    = cloudFromImage(movingim, scalefilter);
    % [pcfix, imgfix]    = cloudFromImage(fixedim, scalefilter);
    
    [optimizer,metric] = imregconfig("multimodal");
    % optimizer.MaximumIterations = 600;
    % metric.NumberOfHistogramBins = 10;
    imfixsharp = imsharpen(fixedim, 'Radius', scalefilter);
    tformout   = imregtform(movingim, fixedim, 'affine', optimizer,metric);
    
    rfixed  = imref2d(size(fixedim)); 
    newim   = imwarp(movingim, rmoving, tformout, 'OutputView',rmoving);
    fimshow = uint8(255 * single(fixedim)/single(quantile(fixedim,0.999,'all')));
    imshowpair(fimshow, newim)
    %%


    isample    = 32;
    istart     = idxsamples(isample) - 5;
    iend       = idxsamples(isample) + 5;
    ismcurr    = pcsamppre(:, 2) > istart & pcsamppre(:, 2) < iend;
    pcsampcurr = pcsamppre(ismcurr, [1 3]);

    iatcurr     = pcatlaswp(:, 2) > istart & pcatlaswp(:, 2) < iend;
    pcatlascurr = pcatlaswp(iatcurr, [1 3]);

    [yreg,bfit] = pcregisterBCPD(pcatlascurr, pcsampcurr, 'TransformType','Affine',...
    'BCPDPath', bcpd, 'OutlierRatio', 0.1, Beta = 15,...
    Gamma = 1, Verbose = false, ConvergenceTolerance=1e-6);
    
    clf; 
    subplot(1,2,2); hold on;
    if size(pcatlascurr,2)==2
        scatter(pcsampcurr(:, 1), pcsampcurr(:,2) ,5, 'k');
        scatter(yreg(:, 1), yreg(:,2), 5, 'r');
        axis equal;
    else
        scatter3(pcsampcurr(:, 1), pcsampcurr(:,2), pcsampcurr(:,3) ,5, 'k');
        scatter3(yreg(:, 1), yreg(:,2), yreg(:,3),5, 'r');
        axis equal;
    end
    subplot(1,2,1); hold on;
    if size(pcatlascurr,2)==2
        scatter(pcsampcurr(:, 1), pcsampcurr(:,2) ,5, 'k');
        scatter(pcatlascurr(:, 1), pcatlascurr(:,2), 5, 'r');
        axis equal;
    else
        scatter3(pcsampcurr(:, 1), pcsampcurr(:,2), pcsampcurr(:,3) ,5, 'k');
        scatter3(pcatlascurr(:, 1), pcatlascurr(:,2), pcatlascurr(:,3),5, 'r');
        axis equal;
    end
    %%

end
%==========================================================================
resall       = nan(Nslices, 1);

slicedist    = 20;

fprintf('Refining slices... ');
msg = []; slicetic = tic;
for islice = 1:Nslices
    %%
    islice = 900;
    isampkeep  = abs(pcsamp(:, iaxuse) - islice) < slicedist;
    Xslicecurr = pcsamp(isampkeep, axkeep);

    % we find the x-y points in the original atlas space
    iatlascurr = abs(pcatlastrans(:,iaxuse) - islice) < slicedist;
    Xatlascurr = pcatlastrans(iatlascurr, axkeep);
    % Xatlascurr = Xatlascurr;



    [yreg,bfit] = pcregisterBCPD(Xatlascurr, Xslicecurr, 'TransformType','Rigid',...
    'BCPDPath', bcpd, 'OutlierRatio', 0.01, ...
    'Gamma', 1, Verbose = false, ConvergenceTolerance=1e-6);

   

    Xslicecurr(:, 2) = randn(size(Xslicecurr, 1), 1) *slicedist/5;
    Xatlascurr(:, 2) = randn(size(Xatlascurr, 1), 1) *slicedist/5;

    pcatlascurr = pointCloud(Xatlascurr);
    Natcurr     = nnz(iatlascurr);
    ndown       = max(6, ceil(Natcurr/1000));
    % pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample', ndown);
    pcatlascurr = pcdownsample(pcatlascurr, 'random', min(1,1500/Natcurr));
    
    pcslicecurr = pointCloud(Xslicecurr);
    Nslicecurr  = pcslicecurr.Count;
    ndown2      = max(6, ceil(Nslicecurr/2000));
    % pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', ndown2);
    pcslicecurr = pcdownsample(pcslicecurr, 'random',  min(1,6000/Nslicecurr));

    meanslice = mean(pcslicecurr.Location(:,[3 1]));
    meanatlas = mean(pcatlascurr.Location(:,[3 1]));
    Tdiff     = -meanatlas + meanslice;


    [R,T,data2] = icp(pcslicecurr.Location(:,[3 1]), pcatlascurr.Location(:,[3 1]) + Tdiff, 100, 10,1, 1e-6);
    res         = mean(min(pdist2(data2', pcslicecurr.Location(:,[3 1])),[],1));
    tformcurr   = rigidtform2d(R, T + Tdiff');
    
    % Options   = struct('Registration','Rigid', 'Verbose', false,'TolP', 1e-4);
    % [~, M]    = ICP_finite_2d(pcslicecurr.Location(:,[3 1]), data2', Options);
    % tformaff = rigidtform2d(M);

    if isaffine
        Options   = struct('Registration', texttrans, 'Verbose', false,'TolP', 1e-4);
        [~, M]    = ICP_finite_2d(pcslicecurr.Location(:,[3 1]), data2', Options);
        tformcurr = affinetform2d(tformcurr.A * M);
    end

    Tmat                = tformcurr.invert;
    if isaffine
        tformslices(islice) = affinetform2d(Tmat.A * tglobal.A);
    else
        tformslices(islice) = rigidtform2d(Tmat.A * tglobal.A);
    end
    resall(islice)      = res;
    %-----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d done. Took %2.2f s. Mean slice error %2.3f. \n', islice, ...
        Nslices, toc(slicetic), median(resall(1:islice)));
    fprintf(msg);
    %-----------------------------------------------------------------------
    % datatest = tformcurr.transformPointsInverse(pcslicecurr.Location(:,[3 1]));
    % % datatest2 = tformaff.transformPointsInverse(pcslicecurr.Location(:,[3 1]));
    % % [tformcod, regcpd] = pcregistercpd(pcatlascurr, pcslicecurr, ...
    % %     'Transform','Affine', 'Tolerance',1e-6,'OutlierRatio',0.0,'Verbose',true,'MaxIterations',100);
    % % regcpd = pctransform(pcslicecurr, tformcod.invert);
    % 
    % clf;
    % subplot(1,3,1)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % subplot(1,3,2)
    % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    %     datatest(:,1), datatest(:,2), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % title('icp')
    % % subplot(1,3,3)
    % % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
    % %     datatest2(:,1), datatest2(:,2), 'k.')
    % % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % % title('cpd')
    % % pause;
    %-----------------------------------------------------------------------
end

end