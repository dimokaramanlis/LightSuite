function tformout = alignCordAtlasToSample(pcatlas, pcsamp, tformslices, trtype)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

bcpd    = 'C:\Users\karamanl\Documents\GitHub\bcpd\win\bcpd.exe';
Ntargetatlas = 1e5;
Ndownatlas   = max(6, floor(size(pcatlas,1)/Ntargetatlas));
pcdownatlas  = pcdownsample(pointCloud(pcatlas), 'nonuniformGrid', Ndownatlas);
pcdownatlas  = pcdownatlas.Location;
pcsamp2      = tranformCordPointsSlices(pcsamp, tformslices);
Ndownsample  = max(6, floor(size(pcsamp,1)/Ntargetatlas));
pcdownsamp   = pcdownsample(pointCloud(pcsamp2), 'nonuniformGrid', Ndownsample);
pcdownsamp    = pcdownsamp.Location;

[yreg,bfit] = pcregisterBCPD(pcdownatlas, pcdownsamp, 'TransformType','Rigid',...
    'BCPDPath', bcpd, 'OutlierRatio', 0.01, ...
    Gamma = 1, Beta = 15, Verbose = false, ConvergenceTolerance=1e-8);
tformout = bfit;
%%
%==========================================================================
clf;
targetlevel = 500;
subplot(1,2,1)
scatter3(pcdownsamp(:,1), pcdownsamp(:,2), pcdownsamp(:,3), 1, 'k');hold on;
scatter3(yreg(:,1), yreg(:,2), yreg(:,3), 1, 'r'); 
axis equal;
subplot(1,2,2);
isubset   = abs(pcdownsamp(:, 3)-targetlevel) < 150;
isubsetat = abs(yreg(:, 3)-targetlevel) < 150;

scatter(pcdownsamp(isubset,1), pcdownsamp(isubset,2), 1, 'k');hold on;
scatter(yreg(isubsetat,1), yreg(isubsetat,2), 1, 'r'); 
axis equal;

%==========================================================================
end