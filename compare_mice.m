
mouseIds   = {'YX002', 'YX004','YX016', 'YX017', 'DK025','DK027',...
    'YX001', 'YX003','YX005','YX007','YX009','DK026', 'DK028',...
    'YX011','YX013', 'DK031',...
    'YX012'};

isolo  = [7:13];
ijoint = [1:6];
iobs   = [14:16];

% DK025, DK027 joint

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
% st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
parcelinfo  = readtable(fullfile(allen_atlas_path,'parcellation_to_parcellation_term_membership.csv'));
[defaultdens, defaultinds] = standardAreaDensities();
%%

bkgthres  = 0.5;
groupinds = unique(parcelinfo.parcellation_index);
isub      = contains(parcelinfo.parcellation_term_set_name, 'substructure');
parcelred = sortrows(parcelinfo(isub,:),'parcellation_index');
Ngroups   = numel(groupinds);

Nmice           = numel(mouseIds);
groupstrs       = lower(parcelred.parcellation_term_name);
groupstrsshort  = lower(parcelred.parcellation_term_acronym);

countsall = nan(Ngroups, 2, Nmice, 'single');
bkgsigall = nan(Ngroups, 2, Nmice, 'single');
diametall = nan(Ngroups, 2, Nmice, 'single');

locsall   = cell(Nmice, 1);
for mouseid = 1:Nmice
    atlasptcoords                    = loadMouseAtlasPoints(mouseIds{mouseid});
    cleanatlaspts                    = sanitizeCellCoords(atlasptcoords, av);
    [areacounts, areavols, idcareas] = groupCellsIntoLeafRegions(cleanatlaspts, av, groupinds);

    % diametall(:, 1, :) = accumarray(idcareas{1}(:,2), cleanatlaspts(idcareas{1}(:,1),4),[],@median);
    % diametall(:, 2, :) = accumarray(idcareas{2}(:,2), cleanatlaspts(idcareas{2}(:,1),4),[],@median);
    % % cellstats1 = accumarray(idcareas{1}(:,2), cleanatlaspts(idcareas{1}(:,1),4),[],@(x) {x});
    % % cellstats2 = accumarray(idcareas{2}(:,2), cleanatlaspts(idcareas{2}(:,1),4),[],@median);

    bgcurr                           = loadMouseBackgroundSignal(mouseIds{mouseid});

    bkgsigall(:, :, mouseid)    = bgcurr;

    % here we discard low-intensity areas from cell counting
    areacounts(bgcurr<bkgthres) = nan;
    countsall(:, :, mouseid)    = areacounts;
    locsall{mouseid} = cleanatlaspts;
    fprintf('Loading cells and background from mouse %d/%d\n', mouseid, Nmice);
end
%%
% plot overall stats
countsplot             = countsall;
bkgplot                = bkgsigall;

%%
f = figure;
f.Units = 'centimeters';
f.Position = [2 2 45 16];
Ncols  = ceil(Nmice/2);
p = panel();
p.pack('v', 2)
for ii = 1:2
    p(ii).pack('h', Ncols)
end
p.de.margin = 1;
txtsize = 10;
p.fontsize = txtsize;
p.margin = [11 1 2 2];
p.de.marginleft = 4;

for imouse = 1:Nmice
    
    irow = floor((imouse-1)/Ncols)+1;
    icol = mod(imouse - 1, Ncols) + 1;
    p(irow, icol).select();

    plot(bkgplot(:,1,imouse), bkgplot(:,2,imouse), '.', [0 30], [0 30], 'MarkerSize', 5)
    axis equal; xlim([-1 40]); ylim([-1 40])
    xticks([0 20 40]);yticks([0 20 40]);
    title(sprintf('%s, bkg', mouseIds{imouse}))
    xlabel('Right hemisphere')
    yticklabels([])
    if icol == 1
        ylabel('Left hemisphere')
        yticklabels([0 20 40])
    end
end
%%
f = figure;
f.Units = 'centimeters';
f.Position = [2 2 45 16];
Ncols  = ceil(Nmice/2);
p = panel();
p.pack('v', 2)
for ii = 1:2
    p(ii).pack('h', Ncols)
end
p.de.margin = 1;
txtsize = 10;
p.fontsize = txtsize;
p.margin = [11 1 2 2];
p.de.marginleft = 4;

for imouse = 1:Nmice
    
    irow = floor((imouse-1)/Ncols)+1;
    icol = mod(imouse - 1, Ncols) + 1;
    p(irow, icol).select();

    plot(log10(countsplot(:,1,imouse)), log10(countsplot(:,2,imouse)), '.', [0 10], [0 10], 'MarkerSize', 5)
    axis equal; xlim([-0.5 6]); ylim([-0.5  6])
    xticks([0 3 6]);yticks([0 3 6]);
    title(sprintf('%s, log10(counts)', mouseIds{imouse}))
    xlabel('Right hemisphere')
    yticklabels([])
    if icol == 1
        ylabel('Left hemisphere')
        yticklabels([0 3 6])
    end
end

%%
countsuse = countsall;
signaluse = bkgsigall;
countsuse(bkgsigall < 0.5) = nan;
signaluse(signaluse < 0.5) = nan;

[fullcounts, fullsignal, fullvols, fullnames, fullinds] = ...
    reorganizeAreas(countsuse, signaluse, areavols, parcelinfo, 'substructure');

[groupcounts, groupsignal, groupvols, groupnames] = ...
    reorganizeAreas(countsuse, signaluse, areavols, parcelinfo, 'structure');

[coarsecounts, coarsesignal, coarsevols] = ...
    reorganizeAreas(countsall, signaluse, areavols, parcelinfo, 'division');

%%
f = figure;
f.Units = 'centimeters';
f.Position = [2 2 45 16];
Ncols  = ceil(Nmice/2);
p = panel();
p.pack('v', 2)
for ii = 1:2
    p(ii).pack('h', Ncols)
end
p.de.margin = 1;
txtsize = 10;
p.fontsize = txtsize;
p.margin = [10 1 4 2];
p.de.marginleft = 8;

for imouse = 1:Nmice
    
    irow = floor((imouse-1)/Ncols)+1;
    icol = mod(imouse - 1, Ncols) + 1;
    p(irow, icol).select();

    indkeep    = ismember(defaultinds, fullinds);
    densplot   = defaultdens(indkeep, :);
    [~, isort] = sort(fullinds, 'ascend');
    currdensities = (fullcounts(isort,imouse)./fullvols(isort));
    currdensities(currdensities<1e-2) = nan;

    plot(log10(densplot(:,1)), log10(currdensities),'.', 'MarkerSize', 4)
    axis square; xlim([4.2 6]); ylim([-1 5])
    xticks([4.5 5 5.5 6]);yticks([0 3 5]);
    title(mouseIds{imouse})
    xlabel('log expected (cells/mm3)')
    yticklabels([])
    if icol == 1
        ylabel('log measured (cells/mm3)')
        yticklabels([0 3 5])
    end
    [rho, pval] = corr(currdensities,...
        densplot,'Type','Spearman', 'Rows','pairwise');
    % [rhoc, pval] = corr((fullcounts(indkeep,imouse)),...
    %     fullsignal(indkeep,imouse)/mean(fullsignal(indkeep,imouse)),'Type','Spearman');
    text(6, 0.2, sprintf('rho_{neu}  = %.2f', rho(1)), 'HorizontalAlignment', 'right')
    text(6, -0.5, sprintf('rho_{exc} = %.2f', rho(2)), 'HorizontalAlignment', 'right')

    glmfit()


end
%%

f = figure;
f.Units = 'centimeters';
f.Position = [2 2 45 16];
Ncols  = ceil(Nmice/2);
p = panel();
p.pack('v', 2)
for ii = 1:2
    p(ii).pack('h', Ncols)
end
p.de.margin = 1;
txtsize = 10;
p.fontsize = txtsize;
p.margin = [10 1 4 2];
p.de.marginleft = 8;

for imouse = 1:Nmice
    
    irow = floor((imouse-1)/Ncols)+1;
    icol = mod(imouse - 1, Ncols) + 1;
    p(irow, icol).select();

    indkeep        = find(fullvols>1e-4);

    plot((fullcounts(indkeep,imouse)./fullvols(indkeep)),...
        fullsignal(indkeep,imouse)/mean(fullsignal(indkeep,imouse)),...
        '.', 'MarkerSize', 5)
    axis square; xlim([-100 5000]); ylim([-0.1 3])
    xticks([0 2500 5000 ]);yticks([0 1 2 3 ]);
    title(mouseIds{imouse})
    xlabel('Area cell density (cells/mm3)')
    yticklabels([])
    if icol == 1
        ylabel('Area bkg signal (norm.)')
        yticklabels([0 1 2 3])
    end
    [rho, pval] = corr((fullcounts(indkeep,imouse)./fullvols(indkeep)),...
        fullsignal(indkeep,imouse)/mean(fullsignal(indkeep,imouse)),'Type','Spearman');
    [rhoc, pval] = corr((fullcounts(indkeep,imouse)),...
        fullsignal(indkeep,imouse)/mean(fullsignal(indkeep,imouse)),'Type','Spearman');
    text(5000, 0.3, sprintf('rho_{dens}  = %.2f', rho), 'HorizontalAlignment', 'right')
    text(5000, 0.08, sprintf('rho_{count} = %.2f', rhoc), 'HorizontalAlignment', 'right')

end


%%

ilayer6a  = find(contains(lower(fullnames(:,1)),'layer 6a'));
ilayer6b  = find(contains(lower(fullnames(:,1)),'layer 6b'));
ilayer5   = find(contains(lower(fullnames(:,1)),'layer 5'));
ilayer23  = find(contains(lower(fullnames(:,1)),'layer 2/3'));
ilayer4   = find(contains(lower(fullnames(:,1)),'layer 4'));
ilayer1   = find(contains(lower(fullnames(:,1)),'layer 1'));
layercell = {ilayer1, ilayer23, ilayer4, ilayer5, ilayer6a, ilayer6b};
layertitles = {'layer 1', 'layer 2/3','layer 4','layer 5','layer 6a','layer 6b'};

f = figure;
f.Units = 'centimeters';
f.Position = [2 2 45 20];
%%
p = panel;
p.pack('h', numel(layercell))
p.de.marginleft = 25;
p.margin = [22 15 2 5];
idslayers = cat(1, layercell{:});
normfacs     = sum(fullcounts(idslayers,:))./sum(fullvols(idslayers));
alldensities = fullcounts./fullvols./normfacs;
% 
% normfacs     = sum(fullsignal(idslayers,:));
% alldensities = fullsignal./normfacs;

maxc         = quantile(alldensities(cat(1, layercell{:}), :), 0.999,'all');
minc         = quantile(alldensities(cat(1, layercell{:}), :), 0.001,'all');

for ii = 1:numel(layercell)
    il = layercell{ii}; 
    p(ii).select();
    [~, isort] = sort(fullnames(il,2));
    currdensities = alldensities(il(isort), :);
    imagesc(currdensities,[minc maxc])
    yticks(1:numel(il))
    yticklabels(fullnames(il(isort),2))
    ax = gca; ax.YDir = 'reverse'; ax.Colormap = magma;
    axis tight;
    title(layertitles{ii})
    xlabel('Mice'); xticks([1 Nmice])
end



%%

groupdensities  = groupcounts./groupvols;
icortex         = contains(lower(groupnames(:,3)),'isocortex');
cortexdensities = sum(groupcounts(icortex, :))./sum(groupvols(icortex,:));

groupdensities = groupdensities./cortexdensities;
% groupdensities = groupsignal./median(groupsignal(icortex,:));

[unnames,~,iun] = unique(groupnames(:,3));

coarserel = coarsecounts./sum(coarsecounts);
f = figure;
f.Units = 'centimeters';
fw = 30;
fh = 22;
f.Position = [1 1 fw fh];
p = panel();
p.pack( 4, 3)

txtsize = 10;
p.margin = [18 16 1 3];
p.fontsize = txtsize;
p.de.margintop = 18;

% maxc = round(maxc/100)*100;
ytext = 'Density (cortex norm.)';

isolo  = [7:13];
ijoint = [1:6];
% isolo  = [7:11];
% ijoint = [1:4];
psurprise = surpriseKruskalWallis(groupdensities(:, ijoint), groupdensities(:, isolo));

for ii = 1:numel(unnames)
    icol = floor((ii-1)/4) + 1;
    irow = mod(ii-1, 4) + 1;

    ishow = iun == ii;

    [acrosort,isort] = sort(groupnames(ishow,2));
    meanjoint = mean(groupdensities(ishow, ijoint), 2);
    meansolo  = mean(groupdensities(ishow, isolo), 2);
    stdall    = std(groupdensities(ishow, :), [], 2);
    jointdens = groupdensities(ishow, ijoint);
    solodens  = groupdensities(ishow, isolo);
    obsdens  = groupdensities(ishow, iobs);

    indexplot  = (meanjoint - meansolo)./stdall;
    indexplot  = psurprise(ishow);

    maxc = ceil(max([solodens jointdens],[],'all'));

    p(irow, icol).select(); cla;

    plot(1:nnz(ishow),jointdens(isort,:), 'r',1:nnz(ishow),solodens(isort,:), 'b')
    plot(1:nnz(ishow), indexplot(isort), 'k')
    line([1 nnz(ishow)], [0 0], 'LineStyle','--')
    % plot(1:nnz(ishow),jointdens(isort,:), 'r',1:nnz(ishow),solodens(isort,:), 'b',...
    %     1:nnz(ishow),obsdens(isort,:), 'g')
    ylabel(ytext)
    xticklabels(acrosort)
    xticks(1:nnz(ishow))
    ylim([0 maxc]); xlim([0.5 nnz(ishow)+1])
    yticks([0 maxc/2 maxc]);
    ylim([-0.11 2.5]); ylabel('k-w surprise');yticks([-2 -1 0 1 2]);
    ax = gca; ax.Box = 'off'; ax.XTickLabelRotation = 90;
    ymax = ylim;ymax = ymax(2);
        text(1, ymax*0.9, unnames{ii})

    if ii == 5
        text(  nnz(ishow), ymax*0.8,sprintf('Solo (N = %d)', numel(isolo)),...
            'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','b')
        text(  nnz(ishow), ymax*0.7,sprintf('Joint (N = %d)', numel(ijoint)),...
            'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','r')
         % text(  nnz(ishow), maxc*0.6,sprintf('Obs. (N = %d)', numel(iobs)),...
         %    'HorizontalAlignment','right', 'FontSize', txtsize-1, 'Color','g')
    end

end
%%
%=========================================================================a
pathsave = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\Presentations\20250324_seminar';
p.export(fullfile(pathsave,'cell_counts.pdf'), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp');

%%
jointmap = atlasCellMap(locsall(ijoint), size(av), 5,  sum(groupcounts(icortex, ijoint)));
solomap  = atlasCellMap(locsall(isolo),  size(av), 5,  sum(groupcounts(icortex, isolo)));
jointout = atlasCellMap(locsall(iobs),   size(av), 5,  sum(groupcounts(icortex, iobs)));
%%
cf = figure('Position',[50 50 1400 800]);
p = panel();
p.pack('h', 3);
p.de.margin = 1;
cmax = 0.01;

dpsave = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\figures\sec_cell_avg';
for islice = 1:5:size(solomap,1);
    soloim = squeeze(solomap(islice,:,:));
    atlasim = single(squeeze(av(islice,:,1:size(av,3)/2)));
    jointim = squeeze(jointmap(islice,:,:));
    jointoutim = squeeze(jointout(islice,:,:));
    
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    p(1).select();cla;
    imagesc(soloim, [0 cmax]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('task solo')
    p(2).select();cla;
    imagesc(jointim, [0 cmax]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('task with other')
    p(3).select();cla;
    imagesc(jointoutim, [0 cmax]); ax = gca; ax.Visible = 'off'; axis equal; axis tight;
    ax.Title.Visible = 'on';
    ax.YDir = 'reverse'; ax.Colormap = flipud(gray);
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',1)
    title('observing other')
    savepngFast(cf, dpsave, sprintf('slice_%03d', islice),300, 2)

end
%%
Nareas     = size(groupcounts, 1);
xcount     = log(sum(groupcounts))';
xother     = zeros(Nmice, 1);
xtask      = zeros(Nmice, 1);
xsex       = zeros(Nmice, 1);

xother([ijoint iobs]) = 1;
xtask([ijoint isolo]) = 1;
xsex(contains(mouseIds, 'DK02')) = 1;

XX     = [xother xsex xcount];

coeffs = nan(Nareas, size(XX, 2));

for iarea = 1:Nareas
    yy     = groupcounts(iarea,:)';
    bcoeff = glmfit(XX, yy, "poisson","LikelihoodPenalty","jeffreys-prior");
    coeffs(iarea, :) = bcoeff(2:end);
end