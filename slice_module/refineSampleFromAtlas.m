function [tformslices, resall] = refineSampleFromAtlas(pcatlas, pcsamp, tformrigid, atlasframe)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

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

xadd = [atlasframe(1) 0 atlasframe(2)]/2;
% Nptsdata       = size(Xlocs, 1);
% Nptsatlas      = pcatlas.Count;
% Nptsperslice   = ceil(10000/Nslices);
resall                  = nan(Nslices, 1);
tformslices(Nslices, 1) = rigidtform2d;
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
    pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample', ndown, 'PreserveStructure',true);

    pcslicecurr = pointCloud(Xslicecurr);
    Nslicecurr  = pcslicecurr.Count;
    ndown2      = max(6, ceil(Nslicecurr/1000));
    pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', ndown2, 'PreserveStructure',true);

    [R,T,data2] = icp(pcslicecurr.Location(:,[3 1]), pcatlascurr.Location(:,[3 1]), 100, 10, 1, 1e-6);
    res         = mean(min(pdist2(data2', pcslicecurr.Location(:,[3 1])),[],1));

    tformcurr   = rigidtform2d(R, T);
    Tmat        = tformcurr.invert;
    tformslices(islice) = rigidtform2d(Tmat.A * tglobal.A);
    resall(islice) = res;
    %-----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d done. Took %2.2f s. Median error %2.3f. \n', islice, ...
        Nslices, toc(slicetic), median(resall(1:islice)));
    fprintf(msg);
    %-----------------------------------------------------------------------
    % datatest = tformcurr.transformPointsInverse(pcslicecurr.Location(:,[3 1]));
    % % [tformcod, regcpd] = pcregistercpd(pcslicecurr, pcatlascurr, ...
    % %     'Transform','Affine', 'Tolerance',1e-6,'OutlierRatio',0.0,'Verbose',true,'MaxIterations',100);
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
    % %     regcpd.Location(:,3), regcpd.Location(:,1), 'k.')
    % % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % % title('cpd')
    % pause;
    %-----------------------------------------------------------------------
end

end