function [cpsample, cpatlas] = triageAndMatchClouds(sampcl, atlascl, transinit, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
%%
Nsample         = sampcl.Count;
Natlas          = atlascl.Count;
rng(1);
%--------------------------------------------------------------------------
Nneighatlas  = ceil(Natlas/5e4);
if Nneighatlas >= 6
    atlascluse = pcdownsample(atlascl,'nonuniformGrid',Nneighatlas,'PreserveStructure',true);
    atlascluse = atlascluse.Location;
else
    atlascluse = atlascl.Location;
end
Nneighsamp = ceil(Nsample/1e4);
if Nneighsamp >= 6
    sampcluse = pcdownsample(sampcl,'nonuniformGrid',Nneighsamp,'PreserveStructure',true);
    sampcluse = sampcluse.Location;
else
    sampcluse = sampcl.Location;
end
sampcluse = transinit.transformPointsForward(sampcluse);
%--------------------------------------------------------------------------
[~, Dclosest]   = knnsearch(atlascluse, sampcluse, 'K', 1);
[~, Dclosest2]  = knnsearch(sampcluse, atlascluse, 'K', 1);
irem  = Dclosest(:,  1) >50;
irem2 = Dclosest2(:, 1) >50;
sampcluse  = sampcluse(~irem, :);
atlascluse = atlascluse(~irem2, :);
% 
% clf; hold on;
% scatter3(atlascluse(:,1), atlascluse(:,2), atlascluse(:,3),1, "filled", 'r')
% scatter3(sampcluse(:,1), sampcluse(:,2), sampcluse(:,3), 1,"filled", 'k')
% axis equal; axis tight; ax = gca; ax.Visible = 'off';
%==========================================================================
% the transform is being fit to find matching pairs
[yreg,bfit] = pcregisterBCPD(atlascluse, sampcluse,'TransformType','AffineNonrigid',...
    'BCPDPath', params.bpcdpath, 'OutlierRatio', 0.01, ...
    Gamma = 0.1, Lambda = 3, Beta = 15, ...
    Verbose = false, ConvergenceTolerance=1e-6, NormalizeCommon = true);

% clf; hold on;
% scatter3(yreg(:,1), yreg(:,2), yreg(:,3),1, "filled", 'r')
% scatter3(sampcluse(:,1), sampcluse(:,2), sampcluse(:,3), 1,"filled", 'k')
% axis equal; axis tight; ax = gca; ax.Visible = 'off';
%==========================================================================
% we finally select points
[Idx, Dnew] = knnsearch(yreg, sampcluse, 'K', 40);
ikeep       = Dnew(:,1) < 20 & median(Dnew,2) < 50;

cpsample   = transinit.transformPointsInverse(sampcluse(ikeep, :));
cpatlas    = atlascluse(Idx(ikeep, 1), :);
%==========================================================================
%%
% clf; 
% subplot(1,2,1); hold on;
% scatter3(cpatlas(:,1), cpatlas(:,2), cpatlas(:,3),1, "filled", 'r')
% scatter3(cpsample(:,1), cpsample(:,2), cpsample(:,3), 1,"filled", 'k')
% axis equal; axis tight; ax = gca; ax.Visible = 'off';
% taff = fitAffineTrans3D(cpatlas, cpsample);
% cpaff = taff.transformPointsForward(cpatlas);
% subplot(1,2,2); hold on;
% scatter3(cpaff(:,1), cpaff(:,2), cpaff(:,3),1, "filled", 'r')
% scatter3(cpsample(:,1), cpsample(:,2), cpsample(:,3), 1,"filled", 'k')
% axis equal; axis tight; ax = gca; ax.Visible = 'off';

%%

end