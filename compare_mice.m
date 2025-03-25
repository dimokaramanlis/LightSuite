mouseIds   = {'DK025', 'DK027', 'YX002', 'YX004','DK026', 'DK028', ...
    'YX001', 'YX003','YX005','YX011','YX012', 'DK031'};
%settingids = {'joint', 'joint', 'solo', 'solo', 'solo', 'joint','observ.'};

isolo = [7 8 9];
ijoint = [3 4];

% isolo  = [5 6 7 8 9];
% ijoint = [1 2 3 4];

iobs = [10 12];

% DK025, DK027 joint

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
%%

bkgthres = 0.5;
Ngroups  = max(av, [], 'all');
Nmice    = numel(mouseIds);
groupstrs       = lower(table2cell(st(1:Ngroups,4)));
groupstrsshort   = table2cell(st(1:Ngroups,5));

countsall = nan(Ngroups, 2, Nmice, 'single');
bkgsigall = nan(Ngroups, 2, Nmice, 'single');
locsall   = cell(Nmice, 1);
for mouseid = 1:Nmice
    atlasptcoords               = loadMouseAtlasPoints(mouseIds{mouseid});

    cleanatlaspts                    = sanitizeCellCoords(atlasptcoords, av);
    [areacounts, areavols, idcareas] = groupCellsIntoLeafRegions(cleanatlaspts, av);
    bgcurr                           = loadMouseBackgroundSignal(mouseIds{mouseid});

    bkgsigall(:, :, mouseid)    = bgcurr;

    % here we discard low-intensity areas from cell counting

    areacounts(bgcurr<bkgthres) = nan;
    countsall(:, :, mouseid)    = areacounts;
    locsall{mouseid} = cleanatlaspts;
    % 
    % avplot = ndSparse.build(cleanatlaspts(:,[2 1 3]), 1,size(av));
    % avplot = single(full(avplot));
    % avgauss = imgaussfilt3(avplot, 5);

end
%%


jointmap = atlasCellMap(locsall(1:4), size(av), 8);
solomap  = atlasCellMap(locsall(5:9), size(av), 8);
jointout = atlasCellMap(locsall([10 12]), size(av), 8);


%%
cf = figure('Position',[50 50 1400 800]);
p = panel();
p.pack('h', 3);
p.de.margin = 1;

dpsave = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\figures\sec_cell_avg';
for islice = 1:5:size(solomap,1);
    soloim = squeeze(solomap(islice,:,:));
    atlasim = single(squeeze(av(islice,:,1:size(av,3)/2)));
    jointim = squeeze(jointmap(islice,:,:));
    jointoutim = squeeze(jointout(islice,:,:));
    
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    p(1).select();cla;
    imagesc(soloim, [0 0.003]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('task solo')
    p(2).select();cla;
    imagesc(jointim, [0 0.003]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('task with other')
    p(3).select();cla;
    imagesc(jointoutim, [0 0.003]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('observing other')
    savepngFast(cf, dpsave, sprintf('slice_%03d', islice),300, 2)

end
%%
countsuse = countsall;
signaluse = bkgsigall;
countsuse(bkgsigall < 0.5) = nan;
signaluse(signaluse < 0.5) = nan;

[newcounts, newsignal, newvols, stnew, namesuse] = reorganizeAreas(countsall, signaluse, areavols, st, 2);

[coarsecounts, coarsesignal, coarsevols, stcoarse] = reorganizeAreas(countsall, signaluse, areavols, st, 1);

newnames = lower(stnew.name);
Nareasfin = size(stnew, 1);


% densitiesall = newsignal;
% densitiesall = densitiesall./(coarsesignal(1,:));
densitiesall  = newcounts./newvols;
%%
[newcounts, newsignal, newvols, stnew, namesuse] = reorganizeAreas(countsall, signaluse, areavols, st, 3);
densitiesall  = newcounts./newvols;

icortex = contains(lower(namesuse(:,3)),'isocortex');
cortextcounts = densitiesall(icortex, :);
cortextcounts = cortextcounts./sum(cortextcounts);

% plot(cortextcounts(:,1:9),'Color',[0 0 0 0.5])
imice = 1:9;
xnames = namesuse(icortex, 2);
[xsorted,isort] = sort(xnames);
imagesc(cortextcounts(isort,:)',[0 0.01])
idxuse = 1:5:size(cortextcounts,1);
xticks(idxuse)
ax = gca; ax.Box = 'off'; ax.FontSize = 15;
xticklabels(xsorted(idxuse))
ylabel('Mouse id'); yticks([1 Nmice])
%%
motordens     = densitiesall(contains(lower(newnames),'primary motor'), :);
% 
% coarsedens   = coarsecounts./coarsevols;
% densitiesall = densitiesall./(coarsedens(1,:));


% 
% % densitiesall = squeeze(nanmean(newsignal, 2));
% % densitiesall = densitiesall./(coarsecounts(1,:));
% 
% summotor     = squeeze(sum(densitiesall(...
%     contains(lower(newnames),'primary motor'), :), 1, 'omitnan'));


% densitiesnorm = densitiesall./sqrt(sum(densitiesall.^2,1));
% [aa, bb] = pca(densitiesnorm');
% maxc = max(abs(bb(:,1)))*1.5;
% 
% subplot(1,2,1)
% plot(bb(:,1), bb(:, 2), 'o'); axis equal;xlim([-1 1]*maxc);ylim([-1 1]*maxc)
% text(bb(:,1),bb(:,2),mouseIds)
% subplot(1,2,2)
% plot(bb(:,1), bb(:, 3), 'o'); axis equal;xlim([-1 1]*maxc);ylim([-1 1]*maxc)
% text(bb(:,1),bb(:,3),mouseIds)

signalsolo = median(densitiesall(:,isolo),2);
signaljoint = median(densitiesall(:,ijoint),2);
signalobs   = median(densitiesall(:,iobs),2);

% [~, isort] = sort(signalsolo, 'descend');
% plot(1:Nareasfin, signalsolo(isort), ...
%     1:Nareasfin, signaljoint(isort),...
%     1:Nareasfin, signalobs(isort))
% 
% corr(densitiesall(:,ijoint))
% cmat = corr(densitiesall, 'Type','Spearman');
% betweengroups = cmat(ijoint, isolo);
% 
% withingroups1 = triu(cmat(ijoint, ijoint),1);
% withingroups2 =  triu(cmat(isolo, isolo),1);
% withingroups  = [withingroups1(withingroups1~=0);...
%     withingroups2(withingroups2~=0)];
% [ cmat(isolo, isolo)]


[~, isort] = sort(std(densitiesall(:, [ijoint isolo iobs]),[],2),'descend');
[~, isort] = sort(std([signalsolo signaljoint signalobs],[],2),'descend');


[unnames,~,iun] = unique(namesuse(:,3));


f = figure;
f.Units = 'centimeters';
fw = 30;
fh = 22;
f.Position = [1 1 fw fh];
p = panel();
p.pack( 4, 2)

txtsize = 10;
p.margin = [18 16 1 3];
p.fontsize = txtsize;
p.de.margintop = 18;

maxc = max(densitiesall(:, [isolo ijoint]),[],'all');
% maxc = round(maxc/100)*100;
if maxc < 10
        ytext = 'Density (norm.)';

else
        ytext = 'Density (cells/mm^3)';

end

for ii = 1:numel(unnames)
    icol = floor((ii-1)/4) + 1;
    irow = mod(ii-1, 4) + 1;

    ishow = iun == ii;

    [acrosort,isort] = sort(stnew.acronym(ishow));
    jointdens = median(densitiesall(ishow, ijoint), 2);
    solodens  = median(densitiesall(ishow, isolo), 2);
    obsdens  = median(densitiesall(ishow, iobs), 2);
    jointdens = densitiesall(ishow, ijoint);
    solodens  = densitiesall(ishow, isolo);
    obsdens  = densitiesall(ishow, iobs);


    p(irow, icol).select(); cla;

    plot(1:nnz(ishow),jointdens(isort,:), 'r',1:nnz(ishow),solodens(isort,:), 'b')
    % plot(1:nnz(ishow),jointdens(isort,:), 'r',1:nnz(ishow),solodens(isort,:), 'b',...
    %     1:nnz(ishow),obsdens(isort,:), 'g')
    ylabel(ytext)
    xticklabels(acrosort)
    xticks(1:nnz(ishow))
    ylim([0 maxc]); xlim([0 nnz(ishow)+1])
    yticks([0 1000 2000])
    ax = gca; ax.Box = 'off'; ax.XTickLabelRotation = 90;
    text(1, maxc*0.9, unnames{ii})
    if ii == 5
        text(  nnz(ishow), maxc*0.8,sprintf('Solo (N = %d)', numel(isolo)),...
            'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','b')
        text(  nnz(ishow), maxc*0.7,sprintf('Joint (N = %d)', numel(ijoint)),...
            'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','r')
         % text(  nnz(ishow), maxc*0.6,sprintf('Obs. (N = %d)', numel(iobs)),...
         %    'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','g')
    end

end

%=========================================================================a
pathsave = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\Presentations\20250324_seminar';
p.export(fullfile(pathsave,'cell_counts.pdf'), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp');
%%
iprimary = contains(namesuse(:,3),'Isocortex') ;
isub      = contains(namesuse(:,3),'Midbrain');
ifrontal = contains(newnames,'prelimbic')|...
    contains(newnames,'infralimbic')|...
    contains(newnames,'orbital')|...
    contains(newnames,'secondary motor')|...
    contains(newnames,'cingulate');


p(3).select();

plot(1:nnz(isub),densitiesall(isub, ijoint), 'r',...
    1:nnz(isub),densitiesall(isub, isolo), 'b')
ylabel('Cell count (relative to MOp)')
xticklabels(stnew.acronym(isub))
xticks(1:nnz(isub))
ylim([0 2])
ax = gca; ax.Box = 'off';
title('Primary areas')
p(2).select();

plot(1:nnz(iprimary),densitiesall(iprimary, ijoint), 'r',...
    1:nnz(iprimary),densitiesall(iprimary, isolo), 'b')
ylabel('Cell count (relative to MOp)')
xticklabels(stnew.acronym(iprimary))
xticks(1:nnz(iprimary))
ylim([0 2])
ax = gca; ax.Box = 'off';
title('Primary areas')
ylim([0 2]);xlim([0 nnz(iprimary) + 1])
ax.XTickLabelRotation = 45;
p(1).select();

plot(1:nnz(ifrontal),densitiesall(ifrontal, ijoint), 'r',...
    1:nnz(ifrontal),densitiesall(ifrontal, isolo), 'b')
ylabel('Cell count (relative to MOp)')
xticklabels(stnew.acronym(ifrontal))
xticks(1:nnz(ifrontal))
ylim([0 2]);xlim([0 nnz(ifrontal) + 1])
text(nnz(ifrontal), 0.4, 'Solo (n = 3)',...
    'HorizontalAlignment','right','Color','b', 'FontSize',12)
text(nnz(ifrontal), 0.35, 'Joint (n = 2)',...
    'HorizontalAlignment','right','Color','r', 'FontSize',12)

ax = gca; ax.Box = 'off';
title('Frontal areas')

%%
toplot = squeeze(nanmax(newcounts, [], 2));
toplot = toplot./sum(toplot,1);
[aa, bb] = pca(toplot');
% [bb,aa] = nnmf(toplot', 3);
maxc = max(abs(bb(:,1)))*1.5;

subplot(1,2,1)
plot(bb(:,1), bb(:, 2), 'o'); axis equal;xlim([-1 1]*maxc);ylim([-1 1]*maxc)
text(bb(:,1),bb(:,2),mouseIds)
subplot(1,2,2)
plot(bb(:,1), bb(:, 3), 'o'); axis equal;xlim([-1 1]*maxc);ylim([-1 1]*maxc)
text(bb(:,1),bb(:,3),mouseIds)
%%
iplot = areavols>0;
imagesc(densplot(iplot,:)',[0 6000])

XX = densplot(iplot,:)';

find(any(isnan(XX),1))

[aa,bb] = pca(densplot(iplot,:)')

%%
figure;
for ii = 1:Nmice
    subplot(2, 6, ii)
    plot(bkgsigall(:,1,ii), bkgsigall(:,2,ii), '.', [0 30], [0 30], 'MarkerSize', 5)
    axis equal; xlim([0 30]); ylim([0 30])
    title(sprintf('%s, signal rel. to bkg', mouseIds{ii}))
    xlabel('Right hemisphere')
    ylabel('Left hemisphere')
end
%%
figure;
for ii = 1:Nmice
    subplot(2, 6, ii)
    plot(log10(countsall(:,1,ii)+1), log10(countsall(:,2,ii)+1), '.', [0 30], [0 30], 'MarkerSize', 5)
    axis equal; xlim([0 6]); ylim([0 6])
    title(sprintf('%s, log10 cell counts', mouseIds{ii}))
    xlabel('Right hemisphere')
    ylabel('Left hemisphere')
end

%%

irem = areavols==0 | strcmp(groupstrs, 'root');
groupstrsuse = groupstrs(~irem, :);
groupstrsshortuse = groupstrsshort(~irem, :);
densitiesuse = countsall(~irem, :, :);
volsuse      = areavols(~irem, :);

summotor     = squeeze(sum(densitiesuse(contains(groupstrsuse,'primary motor'), :, :), [1 2]));
% densitiesuse = squeeze(nanmean(densitiesuse./volsuse, 2))./summotor';
% densitiesuse = squeeze(nanmean(densitiesuse./volsuse, 2));
densitiesuse = squeeze(nanmean(densitiesuse, 2));

 
ivis = contains(groupstrsuse,'primary motor')|...
    contains(groupstrsuse,'primary visual')|...
    contains(groupstrsuse,'prelimbic')|...
    contains(groupstrsuse,'infralimbic')|...
    contains(groupstrsuse,'primary auditory')|...
    contains(groupstrsuse,'secondary motor')|...
    contains(groupstrsuse,'cingulate');

ivis = contains(groupstrsuse,'primary motor')|...
    contains(groupstrsuse,'primary visual')|...
    contains(groupstrsuse,'prelimbic');

% ivis = contains(groupstrsuse,'cingulate')|...
%     contains(groupstrsuse,'primary motor');

bb = barh(groupstrsshortuse(ivis), densitiesuse(ivis,:));
cmapcurr = cbrewer('qual', 'Dark2', Nmice);
% cmapcurr = cmapcurr([1 2 5 6 3 4 7 8], :);

for ii = 1:numel(mouseIds)
    bb(ii).FaceColor = cmapcurr(ii, :);
    bb(ii).BarWidth  = 1;
end
ax = gca; ax.YDir = 'reverse'; ax.Box = 'off';
xlabel('Density (norm. to MOp)')
legend(mouseIds)%%
%%

cmat = corr(densitiesuse);
imagesc(cmat, [0 1]);
ax = gca; ax.Colormap = magma;
ax.XTickLabels = mouseIds;
ax.YTickLabels = mouseIds;
title('Correlation (all areas)')


%%
socidx     = (mean(densitiesuse(:, 1:2), 2) - mean(densitiesuse(:, 3:4), 2))./(mean(densitiesuse(:, 1:2), 2) + mean(densitiesuse(:, 3:4), 2));
[~, isort] = sort(socidx, 'descend');

plot(densitiesuse(isort,:))




%%
%%
stmat                            = atlasTreeToMat(st);
Nmax = numel(areacounts);
%%
areavolsuse     = areavols; % in mm3
cells_per_group = areacounts;
irem = cells_per_group<2;
cells_per_group(irem) = [];
areavolsuse(irem) = [];
stgroups = stmat(1:Nmax, :);

% % move one up
% lastleaf   = findfirst(~isnan(stgroups),2,1,'last');
% prevbranch = lastleaf - 1;
% idxprev    = sub2ind(size(stgroups), find(prevbranch>0), prevbranch(prevbranch>0));
% 
% 
% idxprev - 1

groupstrs(irem) = [];
