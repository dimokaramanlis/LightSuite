function transinit = originalSimilarityTform(samppts, atlaspts, params, Tflip)
%UNTITLED Summary of this function goes here

%==========================================================================
% downsample clouds to a reasonable size
rng(1);
Ndownls      = round(samppts.Count/1e4);
Ndowntv      = round(atlaspts.Count/5e4);

ls_cloud_use = pcdownsample(samppts,'nonuniformGrid', Ndownls,'PreserveStructure',true);
tv_cloud_use = pcdownsample(atlaspts,'nonuniformGrid',Ndowntv,'PreserveStructure',true);
tvpoints = tv_cloud_use.Location;
lspoints = Tflip.transformPointsForward(ls_cloud_use.Location);
%==========================================================================

% Perform basic similarity transform
[yreg,bfit] = pcregisterBCPD(tvpoints, lspoints, 'TransformType','Similarity',...
    'BCPDPath', params.bpcdpath, 'OutlierRatio', 0.01, ...
    'Gamma', 1, Verbose = false, ConvergenceTolerance=1e-8);
% % remove points that seem off
[~, Dclosest]  = knnsearch(yreg, lspoints);
[~, Dclosest2] = knnsearch(lspoints, yreg);
irem = Dclosest >25;
irem2 = Dclosest2 >25;
%==========================================================================
% % repeat transform estimation
[yreg,bfit] = pcregisterBCPD(tvpoints(~irem2,:), lspoints(~irem,:), 'TransformType','Similarity',...
    'BCPDPath', params.bpcdpath, 'OutlierRatio', 0.01, ...
    'Gamma', 1, Verbose = false, ConvergenceTolerance=1e-8);
% [~, Dclosest2] = knnsearch(yreg, lspoints);

tformfit  = bfit.invert;
transinit = affinetform3d(tformfit.A*Tflip.A);
%==========================================================================
end