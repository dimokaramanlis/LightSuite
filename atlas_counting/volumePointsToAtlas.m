function finalpts = volumePointsToAtlas(inputpts, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
% at first, points should be resampled to the registration resolution

% first, point dim order becomes image order
outpts = inputpts(:, [2 1 3]);
% then we permute to atlas
outpts   = outpts(:, trstruct.how_to_perm);
%--------------------------------------------------------------------------
% point coordinates should be moved to match atlas space
scalefac = trstruct.regvolsize(trstruct.how_to_perm)./trstruct.ori_size(trstruct.how_to_perm);
outpts   = outpts.*scalefac;
% %--------------------------------------------------------------------------
% % then, points pass through the initial (rigid) transform
% outpts = trstruct.tform_rigid_samp20um_to_atlas_20um_px.transformPointsForward(outpts);
%--------------------------------------------------------------------------
% then, points pass through the inverse bspline transform
outpts = outpts(:, [2 1 3]);
outpts = transformControlPoints(trstruct.tform_bspline_samp20um_to_atlas_20um_px, (outpts-1)*0.02);
%--------------------------------------------------------------------------
% finally, points pass through the inverse affine transform that was fitted to
% control points

% the inversion is already implemented
finalpts = trstruct.tform_affine_samp20um_to_atlas_20um_px.transformPointsForward(1+outpts/0.02);

%--------------------------------------------------------------------------
% and super finally, points should be resampled to match the 10um
% resolution of the atlas
finalpts = 2*finalpts;
%--------------------------------------------------------------------------
% we remove points outside of the annotation volume
% irem0 = any(finalpts<0, 2);
% iremx = finalpts(:,1) > trstruct.atlassize(2);
% iremy = finalpts(:,2) > trstruct.atlassize(1);
% iremz = finalpts(:,3) > trstruct.atlassize(3);
% irem  = irem0 | iremx | iremy | iremz;
% finalpts(irem, :) = [];
%--------------------------------------------------------------------------
% nrand = min(size(finalpts,1), 1e5);
% iplot = randperm(size(finalpts,1),nrand);
% figure;
% plotBrainGrid; hold on;
% scatter3(finalpts(iplot,2),finalpts(iplot,3),finalpts(iplot,1),2,'filled','MarkerFaceAlpha',0.5)

%--------------------------------------------------------------------------
end