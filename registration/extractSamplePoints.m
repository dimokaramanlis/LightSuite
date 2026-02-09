function ptcloud = extractSamplePoints(voluse, thresuse)
%EXTRACTMATCHINGGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here
rng(1);
gradall    = imgradient3(voluse);
sizevol    = size(voluse);
Nmax        = max(sizevol);
batchsize   = ceil(Nmax/3);
Nbatches    = ceil(sizevol/batchsize);
ptsall      = cell(prod(Nbatches), 1);
idx         = 1;
overallmode = mode(voluse, 'all');



for ibatchy = 1:Nbatches(1)
    istarty = (ibatchy - 1) * batchsize + 1;
    iendy   = min(batchsize * ibatchy, sizevol(1));
    ysamp   = istarty:iendy;
    for ibatchx = 1:Nbatches(2)
        istartx = (ibatchx - 1) * batchsize + 1;
        iendx   = min(batchsize * ibatchx, sizevol(2));
        xsamp   = istartx:iendx;
        for ibatchz = 1:Nbatches(3)
            istartz = (ibatchz - 1) * batchsize + 1;
            iendz   = min(batchsize * ibatchz, sizevol(3));
            zsamp   = istartz:iendz;
            gradcurr = gradall(istarty:iendy, istartx:iendx,istartz:iendz);
            volcurr  = voluse(istarty:iendy, istartx:iendx,istartz:iendz);
            volcurr  = medfilt3(volcurr, [1 1 1]);

            indsthres  = randperm(numel(volcurr), min(1e4, numel(volcurr)));
            thresinit  = quantile(volcurr(indsthres), 0.05, 'all') * 2;
            thresinit  = max(thresinit, overallmode);
            gradcurr(volcurr <thresinit) = 0;
            ipts        = find((gradcurr./volcurr)>thresuse & (volcurr > thresinit));
            [rr,cc,dd]  = ind2sub(size(volcurr), ipts);
            X           = [xsamp(cc)',ysamp(rr)',zsamp(dd)'];
            ptsall{idx} = X;
            idx         = idx  + 1;
            % ptcloud    = pointCloud(X);
            % ptcloud    = pcdownsample(ptcloud,'random',0.1,'PreserveStructure',true);
            % pcshow(ptcloud);
            % pause;
        end
    end
end
ptscat  = cat(1, ptsall{:});
iremx   = ptscat(:, 1) < 2 | ptscat(:, 1) > sizevol(2) - 2;
iremy   = ptscat(:, 2) < 2 | ptscat(:, 2) > sizevol(1) - 2;
iremz   = ptscat(:, 3) < 2 | ptscat(:, 3) > sizevol(3) - 2;
irem    = iremz|iremx|iremy;
ptscat  = pointCloud(ptscat(~irem, :));
ptcloud = pcdownsample(ptscat,'random',0.1,'PreserveStructure',true);
% ptcloud = pcdownsample(ptscat,'nonuniformgrid',50,'PreserveStructure',true);
ptcloud = pcdenoise(ptcloud);


% indsthres  = randperm(numel(voluse), min(1e4, numel(voluse)));
% thresinit  = quantile(voluse(indsthres), 0.05, 'all') * 2;
% testout(voluse <thresinit) = 0;
% ipts       = find(testout>thresuse & voluse > thresinit);
% [rr,cc,dd] = ind2sub(size(voluse), ipts);
% X          = [cc,rr,dd];
% ptcloud    = pointCloud(X);
% Nfit       = 5e4;
% pdown      = min(1, Nfit/ptcloud.Count);
% % ptcloud    = pcdownsample(ptcloud,'nonuniformGridSample',pdown,'PreserveStructure',true);
% ptcloud    = pcdownsample(ptcloud,'random',pdown,'PreserveStructure',true);

 % scatter3(ptcloud.Location(:,1),ptcloud.Location(:,2),ptcloud.Location(:,3),2)
 %%
% islice = 200;
% subplot(1,3,1)
% maxc = quantile(voluse, 0.999,'all');
% imtoplot = squeeze(voluse(:,islice,:));
% imagesc(imtoplot, [0 maxc])
% axis image off; ax =gca; ax.Colormap = flipud(gray);
% subplot(1,3,3)
% iplot = abs(ptcloud.Location(:,1) - islice)<15;
% scatter(ptcloud.Location(iplot,3), ptcloud.Location(iplot,2), 8, 'filled','MarkerFaceColor','k')
% ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
% axis equal; xlim([0.5 size(imtoplot, 2)]);
% ylim([0.5 size(imtoplot, 1)])

%%
end

