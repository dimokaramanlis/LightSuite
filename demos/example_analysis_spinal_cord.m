% example analysis scipt for quantifying intensities in spinal cord data
%==========================================================================
% we load the atlas
[tv, av, parcelinfo, segments] = loadSpinalCordAtlas();
% we find unique ids in the annotation volume
[avinds, ~, avic] = unique(av);
Nsegments       = size(segments, 1);
Ngroups         = numel(avinds);
edgebins        = [segments.Start(1:end); segments.End(end)];
segcents        = edgebins(1:end-1) + diff(edgebins)/2;
[~, ~, lenval]  = ind2sub(size(av), (1:numel(av))');
[~, ~, leninds] = histcounts(lenval, edgebins);
%==========================================================================
% we get the volumes of all areas
volumeall       = accumarray([avic leninds], 1, [numel(avinds), Nsegments], @sum);
%==========================================================================
%%
% for all samples, we group signals to the finest level
ichanuse        = 1:2;
% here you have to set up the paths of your samples
dpspinesample   = {'D:\spine_registration\sample1', 'D:\spine_registration\sample2'};
Nchannels       = numel(ichanuse);
Nmice           = numel(dpspinesample);
medianoverareas = nan(Ngroups-1, Nsegments, numel(ichanuse), numel(dpspinesample), 'single');

% volumeall        = zeros([size(av)  numel(dpspinesample)], 'single');
for isample  = 1:Nmice
    volload  = fullfile(dpspinesample{isample}, "lightsuite", "volume_registered");
    imgStack = loadLargeSliceVolume(volload, ichanuse);

    for ichan = 1:numel(ichanuse)
        currvals    = squeeze(imgStack(:,:, ichan,:));
        currvals    = imresize3(currvals, size(av));
        valsgrouped = single(accumarray([avic leninds], currvals(:), [Ngroups, Nsegments], @median));
        relsignal   = (valsgrouped - valsgrouped(1, :))./ valsgrouped(1, :);
        medianoverareas(:, :, ichan, isample) = relsignal(2:end, :);
    end
    fprintf('Sample %d analyzed\n', isample)
end
%==========================================================================
%%
% we then calculate normalization factors with nnmf for each channel
volwts = volumeall(2:end,:)./sum(volumeall(2:end,:),2);
projstrength = nan(numel(ichanuse), numel(dpspinesample));
for ii = 1:numel(ichanuse)
    currmat  = squeeze(medianoverareas(:, :, ii, :));
    norminds = all((currmat > 0) & ~isnan(currmat) & ~isinf(currmat), 3);
    Nuse     = nnz(norminds);
    datacurr = reshape(currmat(repmat(norminds, [1 1 Nmice])), [Nuse Nmice]);
    [projstrength(ii,:), ~] = nnmf(datacurr', 1, "replicates",10);

    structres(ii, 1) = reorganizeSpinalCordAreas(squeeze(medianoverareas(:, :, ii, :)), volumeall(2:end, :), ...
        parcelinfo, avinds(2:end), 'structure');
    divres(ii, 1) = reorganizeSpinalCordAreas(squeeze(medianoverareas(:, :, ii, :)), volumeall(2:end, :), ...
        parcelinfo, avinds(2:end), 'division');
end

projstrength  = projstrength./median(projstrength, 2);
projstrength  = reshape(projstrength, [1 size(projstrength)]);
%==========================================================================
% we then plot our results across structures and divisions
channelNames  = {'DAPI', 'eGFP - Myelin'};
yareas        = structres(1).names(:,2);
yareas        = strrep(yareas, '_', ' ');
xareas        = segments.Segment;
indsx         = 1:2:Nsegments;
%%
% create a nice figure
ffig = figure('Units','centimeters');
fh   = 20; fw   = 35;
ffig.Position = [1 1 fw fh];

p = panel();
p.pack('v', 2);
for ii = 1:2
    p(ii).pack('h', Nchannels)
end
p.de.margintop = 25;
p.de.marginleft = 25;
p.margin = [25 20 2 5];


for ichan = 1:2
    imtoplot = median(structres(ichan).signal./projstrength(:, ichan, :),3);

    imtoplot(isinf(imtoplot)) = 0;
    p(1,ichan).select();
    imagesc(imtoplot)
    pbaspect([2 1 1]); axis tight;
    title(channelNames{ichan})
    yticks(1:size(imtoplot,1));
    yticklabels(yareas)
    xticks(indsx)
    xticklabels(xareas(indsx))
    ax = gca; ax.Colormap = hot;
    ax.XTickLabelRotation = 0; ax.YDir = 'reverse';

    imtoplot = median(divres(ichan).signal./projstrength(:, ichan, :),3);
    divnames = divres(ichan).names(:, 1);

    p(2,ichan).select();
    plot(segcents*20*1e-3, imtoplot')
    pbaspect([2 1 1]); xlim([0 30])
    title(channelNames{ichan})
    xticks(segcents(indsx)*20*1e-3)
    xticklabels(xareas(indsx))
    ylabel('Relative signal strength')
    legend(divnames)
    xlabel('Rostrocaudal position')
end

%%
pathsave = 'D:\';
p.export(fullfile(pathsave,'cord_stats.pdf'), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp');

