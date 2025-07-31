function ptcloud = extractVolumePoints(voluse, sizefilt)
%EXTRACTMATCHINGGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

testout    = stdfilt(voluse, ones(sizefilt*[1,1,1])); 
ptthres    = quantile(testout(testout>1e-3), 0.8, 'all');
Nfit       = 2e4;
ipts       = find(testout>ptthres);
% ipts       = ipts(randperm(numel(ipts), Nfit));
[rr,cc,dd] = ind2sub(size(voluse), ipts);
X          = [cc,rr,dd];
ptcloud    = pointCloud(X);
pdown      = min(1, Nfit/ptcloud.Count);
ptcloud    = pcdownsample(ptcloud,'random', pdown);

 % scatter3(ptcloud.Location(:,1),ptcloud.Location(:,2),ptcloud.Location(:,3),2)
end

