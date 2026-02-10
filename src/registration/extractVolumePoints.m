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
 %%
% islice = 200;
% subplot(1,3,1)
% maxc = quantile(voluse, 0.999,'all');
% imtoplot = squeeze(voluse(:,islice,:));
% imagesc(imtoplot, [0 maxc])
% axis image off; ax =gca; ax.Colormap = flipud(gray);
% subplot(1,3,2)
% imagesc(squeeze(testout(:,200,:)), [0 2*ptthres])
% axis image off; ax =gca; ax.Colormap = flipud(gray);
% subplot(1,3,3)
% iplot = abs(ptcloud.Location(:,1) - islice)<10;
% scatter(ptcloud.Location(iplot,3), ptcloud.Location(iplot,2), 8, 'filled','MarkerFaceColor','k')
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% axis equal; xlim([0.5 size(imtoplot, 2)]);
% ylim([0.5 size(imtoplot, 1)])

%%
end

