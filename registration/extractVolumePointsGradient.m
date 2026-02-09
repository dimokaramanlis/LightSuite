function ptcloud = extractVolumePointsGradient(voluse, sigmause, thresuse)
%EXTRACTMATCHINGGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here
rng(1);

testout    = imgradient3(voluse)./imgaussfilt3(voluse, sigmause,"Padding","symmetric");
indsthres  = randperm(numel(voluse), min(1e4, numel(voluse)));
thresinit  = quantile(voluse(indsthres), 0.05, 'all') * 2;
testout(voluse <thresinit) = 0;
ipts       = find(testout>thresuse & voluse > thresinit);
[rr,cc,dd] = ind2sub(size(voluse), ipts);
X          = [cc,rr,dd];
ptcloud    = pointCloud(X);
% Nfit       = 5e4;
% pdown      = min(1, Nfit/ptcloud.Count);
% % ptcloud    = pcdownsample(ptcloud,'nonuniformGridSample',pdown,'PreserveStructure',true);
% ptcloud    = pcdownsample(ptcloud,'random',pdown,'PreserveStructure',true);

 % scatter3(ptcloud.Location(:,1),ptcloud.Location(:,2),ptcloud.Location(:,3),2)
 %%
% islice = 100;
% subplot(1,3,1)
% maxc = quantile(voluse, 0.999,'all');
% imtoplot = squeeze(voluse(:,islice,:));
% imagesc(imtoplot, [0 maxc])
% axis image off; ax =gca; ax.Colormap = flipud(gray);
% subplot(1,3,2)
% imagesc(squeeze(testout(:,islice,:)), [0 2*ptthres])
% axis image off; ax =gca; ax.Colormap = flipud(gray);
% subplot(1,3,3)
% iplot = abs(ptcloud.Location(:,1) - islice)<15;
% scatter(ptcloud.Location(iplot,3), ptcloud.Location(iplot,2), 8, 'filled','MarkerFaceColor','k')
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% axis equal; xlim([0.5 size(imtoplot, 2)]);
% ylim([0.5 size(imtoplot, 1)])

%%
end

