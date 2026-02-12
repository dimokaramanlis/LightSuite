% this is an example analysis scipt for quantifying intensities in spinal
% cord data


sizeside            = 20; %um
volume_mm3          = (sizeside*1e-3)^3;
dpspineatlas        = 'D:\AllenAtlas\extra_spine\allen_cord_20um_v1.1';
[tv, av, parcelinfo, segments] = loadSpinalCordAtlas(dpspineatlas, 20);
edgebins            = [segments.Start(1:end); segments.End(end)];
[areaidx, ~, icun]  = unique(av);
Ngroups             = numel(areaidx);
volumeoverareas     = single(accumarray(icun,        1, [Ngroups 1], @sum)) * volume_mm3;

[~, ~, lenval]  = ind2sub(size(av), (1:numel(av))');
[~, ~, leninds] = histcounts(lenval, edgebins);
lencents        = edgebins(1:end-1) + diff(edgebins)/2;
Nsegments       = numel(edgebins) - 1;
% iwmall  = contains(parcelinfo.structure_id_path, wmcode);

iwm         = contains(lower(parcelinfo.name), 'white matter');
wmcode      = parcelinfo.structure_id_path{iwm};
igm         = contains(lower(parcelinfo.name), 'gray matter');
gmcode      = parcelinfo.structure_id_path{igm};

wmids       = parcelinfo.id(contains(parcelinfo.structure_id_path, wmcode));
gmids       = parcelinfo.id(contains(parcelinfo.structure_id_path, gmcode));

icurrwm = reshape(ismember(av, wmids), [], 1);
icurrgm = reshape(ismember(av, gmids), [], 1);

atlasgmoverlen = accumarray(leninds(icurrgm), 1, [max(leninds) 1], @sum);
atlaswmoverlen = accumarray(leninds(icurrwm), 1, [max(leninds) 1], @sum);
ratioverlen    = atlasgmoverlen./(atlasgmoverlen +atlaswmoverlen);


%%
ichanuse        = 1:2;
dpspinesample    = {'D:\spine_registration\sample1', 'D:\spine_registration\sample2'};

medianoverareas   = nan(Ngroups, numel(ichanuse), numel(dpspinesample), 'single');
graymattersignal  = nan(Nsegments, numel(ichanuse), numel(dpspinesample), 'single');
whitemattersignal = nan(Nsegments, numel(ichanuse), numel(dpspinesample), 'single');

% volumeall        = zeros([size(av)  numel(dpspinesample)], 'single');
for isample = 1:numel(dpspinesample)
    volload  = fullfile(dpspinesample{isample}, "lightsuite", "volume_registered");
    imgStack  = loadLargeSliceVolume(volload, ichanuse);

    for ichan = 1:numel(ichanuse)
        sidevals  = reshape(imgStack(:,:, ichan,:), [], 1);
        medianoverareas(:, ichan, isample)   = single(accumarray(icun, sidevals, [Ngroups 1], @median));
        graymattersignal(:, ichan, isample)  =...
            accumarray(leninds(icurrgm), sidevals(icurrgm), [Nsegments 1], @median);
        whitemattersignal(:, ichan, isample) = ...
            accumarray(leninds(icurrwm), sidevals(icurrwm), [Nsegments 1], @median);
    end
    backsignal = medianoverareas(1, :, isample);
    medianoverareas(:, :, isample)   = (medianoverareas(:, :, isample) - backsignal)./backsignal;
    graymattersignal(:, :, isample)  = (graymattersignal(:, :, isample) - backsignal)./backsignal;
    whitemattersignal(:, :, isample) = (whitemattersignal(:, :, isample) - backsignal)./backsignal;
end
projstrength = nan(numel(ichanuse), numel(dpspinesample));
for ii = 1:numel(ichanuse)
    [projstrength(ii,:), ~] = nnmf(squeeze(medianoverareas(2:end, ii, :))', 1, "replicates",10);
end
projstrength  = projstrength./median(projstrength, 2);
projstrength  = reshape(projstrength, [1 size(projstrength)]);
%%
toplot = medianoverareas./projstrength;
toplotwm = whitemattersignal./projstrength;
toplotgm = graymattersignal./projstrength;
toplotgm(toplotgm==0) = nan;
toplotwm(toplotwm==0) = nan;

colsgraywhite = [0 0 0; 0.5 0.5 0.5];
channelNames  = {'DAPI', 'eGFP - Myelin', ' '};
% pbnum = 5;
xcord  = (1:size(tv, 3))*sizeside*1e-3;
ycord  = (1:size(tv, 2))*sizeside*1e-3;

xmax = floor(max(xcord)/2)*2;

cf = figure;
cf.Units = "centimeters";
fh = 15; fw = 24;
cf.Position = [1 1 fw fh];

p = panel();
p.pack('v', {0.3 0.35 0.35});
p(2).margintop = 1;
p(3).margintop = 20;

p.margin =[20 20 3 5];


p(1).select(); cla;
imtoshow = squeeze(max(imgStack(:,:,2,:),[],1));
% imtoshow = squeeze(imgStack(60,:,2,:));
imlimscurr = getImageLimits(imgStack,0.01);
imagesc(xcord, ycord, imtoshow, imlimscurr);
xline(edgebins*sizeside*1e-3,'--','Color',[1 0.5 0.5])
for ii = 1:2:size(segments,1)
    text(lencents(ii)*sizeside*1e-3, -0.5, segments.Segment{ii}, ...
        'HorizontalAlignment','center','Color',[1 0.5 0.5])
end

ax = gca; ax.Colormap = gray; ax.Visible = 'off';
ax.Title.Visible = 'on'; 
title('Registered spinal cord sample')
% imagesc(xcord, [], squeeze(imgStack(60,:,2,:)));
axis equal;
% pbaspect([pbnum 1 1]);
xlim([-0.5 xmax]);  ylim([0 max(ycord)]);
for ichan = 1:numel(ichanuse)
    currsig = cat(3, squeeze(toplotgm(:, ichan, :)), squeeze(toplotwm(:, ichan, :)));
    ymax = ceil(max(currsig,[], 'all')/2)*2;
    ymin = 0;
    p(1+ichan).select();cla;
    xline(edgebins*sizeside*1e-3,'--','Color',[1 0.5 0.5])
    line(lencents*sizeside*1e-3, squeeze(toplotgm(:, ichan, :)), 'Color', colsgraywhite(1,:))
    line(lencents*sizeside*1e-3, squeeze(toplotwm(:, ichan, :)), 'Color', colsgraywhite(2,:))
    line(lencents*sizeside*1e-3, median(toplotgm(:, ichan, :), 3), ...
        'Color', colsgraywhite(1,:), 'LineWidth',2)
    line(lencents*sizeside*1e-3, median(toplotwm(:, ichan, :), 3), ...
        'Color', colsgraywhite(2,:), 'LineWidth',2)
    % pbaspect([pbnum 1 1]); 
    ylim([ymin ymax])
    yticks([ymin mean([ymin,ymax]) ymax])
    ylabel('Signal-to-background ratio')
    xlim([-0.5 xmax]); 
    xticks([0 xmax/2 xmax]); xticklabels([])
    title(channelNames{ichan});
    
    text(0, ymin + (ymax-ymin)*0.2, 'Gray matter', 'Color', colsgraywhite(1,:), 'HorizontalAlignment','left')
    text(0, ymin + (ymax-ymin)*0.1, 'White matter', 'Color', colsgraywhite(2,:), 'HorizontalAlignment','left')
    if ichan == 2
        xticks([0 xmax/2 xmax])
        xticklabels([0 xmax/2 xmax])
        xlabel('Rostrocaudal position (mm)')
    end
end


pathsave = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\figures\20251217_spinalcord';
p.export(fullfile(pathsave,'cord_stats.pdf'), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp');


%%
volumeavg              = median(volumeall.*reshape(1./projstrength, [1,1,1,2]), 4);

atlasTreeToMat(parcelinfo);
%%
