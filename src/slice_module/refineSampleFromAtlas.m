function [tformslices, resall] = refineSampleFromAtlas(pcatlas, pcsamp, tformrigid, transformtype)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

rng(1);

pcatlastrans = pctransform(pcatlas, tformrigid);
Natlas       = pcatlastrans.Count;
ikeep        = randperm(Natlas, 1e3);
oldpts       = pcatlas.Location(ikeep, [3 1]);
newpts       = pcatlastrans.Location(ikeep, [3 1]);
tglobal      = fitgeotform2d(newpts, oldpts, 'similarity');
tglobal      = rigidtform2d(tglobal.R, tglobal.Translation);

Xlocs          = pcsamp.Location;
[apun, ~, ic]  = unique(pcsamp.Location(:, 2));

slicedist = median(diff(apun));
Nslices        = numel(apun);

resall                  = nan(Nslices, 1);
switch lower(transformtype)
    case 'rigid'
        tformslices(Nslices, 1) = rigidtform2d;
        isaffine = false;
    case 'affine'
        tformslices(Nslices, 1) = affinetform2d;
        isaffine = true;
        texttrans = 'Affine';
    case 'similarity'
        tformslices(Nslices, 1) = affinetform2d;
        isaffine = true;
        texttrans = 'Size';
end

fprintf('Refining slices... ');
msg = []; slicetic = tic;
for islice = 1:Nslices

    % we find the x-y points in the original atlas space
    Xslicecurr = Xlocs(ic == islice, :);
    apcurr     = apun(islice);
    iatlascurr = abs(pcatlastrans.Location(:,2) - apcurr) < slicedist/2;
    Xatlascurr = pcatlastrans.Location(iatlascurr, :);
    % Xatlascurr = Xatlascurr;

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