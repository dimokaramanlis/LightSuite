function tformout = fitPointCloudsAffine(movingim, fixedim, pxsize, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% pxsize in um
if nargin < 4
    tformtype = 'affine';
else
    tformtype = 'rigid';
end

rng(1);

scalefilter = 100/pxsize;
% [pcmov, imgmov]    = cloudFromImage(movingim, scalefilter);
% [pcfix, imgfix]    = cloudFromImage(fixedim, scalefilter);

[optimizer,metric] = imregconfig("multimodal");
% optimizer.MaximumIterations = 600;
% metric.NumberOfHistogramBins = 10;
imfixsharp = imsharpen(fixedim, 'Radius', scalefilter);
tformout   = imregtform(movingim, imfixsharp, tformtype, optimizer,metric);

% 
% rfixed  = imref2d(size(fixedim)); 
% newim   = imwarp(movingim, rfixed, tformout, 'OutputView',rfixed);
% fimshow = uint8(255 * single(fixedim)/single(quantile(fixedim,0.99,'all')));
% imshowpair(fimshow, newim)


% imtest  = imregister(movingim, imfixsharp, "affine", optimizer,metric);
% tformout = imregtform(imgmov, imgfix, "affine", optimizer,metric);
% 
% subplot(1,2,2)
% imshowpair(fimshow,newim)



%==========================================================================
% using control points - inferior
% 
% [R,T,data2] = icp(pcfix, pcmov, 100, 10, 1, 1e-6);
% pcrigid = data2';
% 
% Options   = struct('Registration','Affine', 'Verbose', true,'TolP', 1e-4, 'TolX', 1e-4);
% [~, M]    = ICP_finite_2d(pcfix, pcmov, Options);
% tformaff = affinetform2d(M);
% pcaffine = tformaff.transformPointsForward(pcmov);
% 
% subplot(1,2,1);
% plot(pcfix(:,1), pcfix(:,2), 'k.', pcmov(:,1), pcmov(:,2), 'r.')
% ax = gca; ax.YDir = 'reverse'; axis equal;
% subplot(1,2,2);
% plot(pcfix(:,1), pcfix(:,2), 'k.', pcaffine(:,1), pcaffine(:,2), 'r.')
% ax = gca; ax.YDir = 'reverse'; axis equal;
% 
% res         = mean(min(pdist2(data2', pcslicecurr.Location(:,[3 1])),[],1));
% tformcurr   = rigidtform2d(R, T);
% 
% 
% 
% 
% movhigh     = spatial_bandpass(single(movingim), scalefilter, 3, 3, false);
% fixhigh     = spatial_bandpass(single(fixedim), scalefilter, 3, 3, false);
% 
% thresfix     = quantile(fixhigh,0.99,'all')/2;
% thresmov     = quantile(movhigh,0.99,'all')/2;
% idxfix          = find(fixhigh>thresfix);
% [rowfix,colfix] = ind2sub(size(fixhigh), idxfix);
% 
% idxmov          = find(movhigh>thresmov);
% [rowmov,colmov] = ind2sub(size(movhigh), idxmov);
% 
% pcmov = pointCloud([colmov rowmov]);
% pcregistercpd()
% 
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
%==========================================================================
end


function [cloudout, imfilt] = cloudFromImage(inputim, scalefilter)

imfilt          = spatial_bandpass(single(inputim), scalefilter, 3, 3, false);
thresfix        = quantile(imfilt,0.99,'all')/2;
idxfix          = find(imfilt>thresfix);
[rowfix,colfix] = ind2sub(size(imfilt), idxfix);
cloudout        = [colfix zeros(size(colfix)) rowfix];
pcout           = pcdenoise(pointCloud(cloudout), "PreserveStructure",true, "NumNeighbors",10);
pcout           = pcdownsample(pcout, "random",0.5, 'PreserveStructure', true);
cloudout        = pcout.Location(:, [1 3]);
end