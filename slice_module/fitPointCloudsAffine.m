function tformout = fitPointCloudsAffine(movingim, fixedim, pxsize)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% pxsize in um

scalefilter = 100/pxsize;
pcmov       = cloudFromImage(movingim, scalefilter);
pcfix       = cloudFromImage(fixedim, scalefilter);


movhigh     = spatial_bandpass(single(movingim), scalefilter, 3, 3, false);
fixhigh     = spatial_bandpass(single(fixedim), scalefilter, 3, 3, false);

thresfix     = quantile(fixhigh,0.99,'all')/2;
thresmov     = quantile(movhigh,0.99,'all')/2;
idxfix          = find(fixhigh>thresfix);
[rowfix,colfix] = ind2sub(size(fixhigh), idxfix);

idxmov          = find(movhigh>thresmov);
[rowmov,colmov] = ind2sub(size(movhigh), idxmov);

pcmov = pointCloud([colmov rowmov]);
pcregistercpd()

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end


function cloudout = cloudFromImage(inputim, scalefilter)

end