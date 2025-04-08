newdataloc   = 'F:\test_manuscript';
cdirfind     = dir(fullfile(newdataloc, '**', 'resfull.mat'));
Nmice        = numel(cdirfind);
%%
iplot      = 2;
datepoints = nan(Nmice, 2);
drpoints   = nan(Nmice, 2);

Ncpoints   = nan(Nmice, 2);



cellsnr     = nan(Nmice, 2);
cellsnrerr  = nan(Nmice, 2);
cellfluo    = nan(Nmice, 2);
cellfluoerr = nan(Nmice, 2);

origplot   = nan;
newplot    = nan;


for ii = 1:Nmice
    % load data
    datamouse = load(fullfile(cdirfind(ii).folder, cdirfind(ii).name));

    % load images
    origvol   = tiffreadVolume(datamouse.oldsignalvol);
    newvol    = tiffreadVolume(datamouse.newsignalvol);
    if ii == iplot
        origplot = origvol;
        newplot  = newvol;
    end


    datevec   = [datamouse.olddate datamouse.newdate];
    datevec   = datevec - datevec(1);
    datevec   = days(datevec);
    datepoints(ii, :) = datevec;

    oldcellxyz = round(datamouse.oldcelllocations(:, 1:3)./[1.44 1.44 5]);
    newcellxyz = round(datamouse.newcelllocs(:, 1:3));
    Ncellsold  = size(oldcellxyz, 1);
    Ncellsnew  = size(newcellxyz, 1);

    Ncpoints(ii, :) = [Ncellsold Ncellsnew];


    indsold   = sub2ind(size(origvol), oldcellxyz(:,2), oldcellxyz(:,1), oldcellxyz(:,3));
    oldsignal = single(origvol(indsold));
    indsnew   = sub2ind(size(newvol), newcellxyz(:,2), newcellxyz(:,1), newcellxyz(:,3));
    newsignal = single(newvol(indsnew));

    cellfluo(ii, :)    = [median(oldsignal) median(newsignal)];
    cellfluoerr(ii, :) = [range(getRobustConfidenceIntervals(oldsignal, 0.95))...
        range(getRobustConfidenceIntervals(newsignal, 0.95))]/2;

    oldf = median(datamouse.oldcelllocations(:, 5));
    newf = median(datamouse.newcelllocs(:, 5));
    oldsnrerr  = getRobustConfidenceIntervals(datamouse.oldcelllocations(:, 5), 0.95);
    newsnarerr = getRobustConfidenceIntervals(datamouse.newcelllocs(:, 5), 0.95);
    cellsnr(ii, :)    = [oldf newf];
    cellsnrerr(ii, :) = [range(oldsnrerr)/2 range(newsnarerr)/2];



    drcurr = [range(datamouse.oldhistogramvals) range(datamouse.newhistogramvals)];

    drpoints(ii, :) = drcurr;

    fprintf('%d/%d mice loaded\n', ii, Nmice)
end
%%
xpix      = size(newplot, 2);

fw = 18;
fh = 14;
ccol = cbrewer('qual','Paired',14);
ccol = ccol([5 6], :);

f = figure;
f.Units = 'centimeters';
f.Position = [2 2 fw fh];

datarange = 250:270;
collims   = [0 3500];
%%
p = panel();
p.pack('v', {0.7 0.3});
p(1).pack('h', 2);
p(2).pack('h', 2)
p(2,1).pack('h', 2)
p(2,2).pack('h', 2)

txtsize = 10;
p.fontsize  = txtsize;
p.de.margin = 1;
p(1).de.marginleft = 2;
p(2).margintop = 1;
p(2,2).marginleft = 16;
p(2,1).de.marginleft = 9;
p(2,2).de.marginleft = 9;

p.margin = [12 7 3 0];

p(1,1).select(); cla;
imagesc(max(origplot(:, :, datarange), [], 3), collims)
text(xpix*0.95, xpix*0.92, '500 um', 'HorizontalAlignment','right','Color','w')
line(xpix*0.95 - [0 500]/1.44, [1 1] *xpix*0.95, 'Color', 'w', 'LineWidth', 2)
axis equal; axis tight; ax = gca; ax.Colormap = copper;
ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off'; 
ax.Title.Visible = 'on'; title('First scan')

p(1,2).select(); cla;
imagesc(max(newplot(:, :, datarange), [], 3), collims)
axis equal; axis tight; ax = gca; ax.Colormap = copper;
ax = gca; ax.YDir = 'reverse'; ax.Visible = 'off';
text(xpix*0.95, xpix*0.92, '500 um', 'HorizontalAlignment','right','Color','w')
line(xpix*0.95 - [0 500]/1.44, [1 1] *xpix*0.95, 'Color', 'w', 'LineWidth', 2)
ax.Title.Visible = 'on'; title(sprintf('Second scan (after %d days)', round(datepoints(iplot,2))))

p(2,1,1).select();cla;
line([0 1e4], [0 1e4], 'Color','k','LineStyle','--')
for ii = 1:Nmice
    currcol = (datepoints(ii,2)>50) + 1;
    line(log10(drpoints(ii,1)),log10(drpoints(ii,2)),'Marker','o', 'Color', 'k',...
        'MarkerFaceColor', ccol(currcol, :))
end
axis equal; ylim([2.95 4]); xlim([2.95 4]);
yticks([3 3.5 4]); xticks([3 3.5 4]);
xlabel('First scan');
ylabel('Second scan')
title('log (dynamic range)')
text(3, 3.9, '~170 days', 'Color',ccol(2,:))
text(3, 3.75, ' ~50 days', 'Color',ccol(1,:))


p(2,1,2).select();cla;
line([0 1e4], [0 1e4], 'Color','k','LineStyle','--')
for ii = 1:Nmice
    currcol = (datepoints(ii,2)>50) + 1;
    line(log10(cellfluo(ii,1)), log10(cellfluo(ii,2)),'Marker','o', 'Color', 'k',...
        'MarkerFaceColor', ccol(currcol, :))
end
axis equal; ylim([2.95 4]); xlim([2.95 4]);
yticks([3 3.5 4]); xticks([3 3.5 4]); yticklabels([])
xlabel('First scan');
title('log (cell fluorescence)')


p(2,2,1).select();cla;
line([0 1e4], [0 1e4], 'Color','k','LineStyle','--')
for ii = 1:Nmice
    currcol = (datepoints(ii,2)>50) + 1;
    line(log10(Ncpoints(ii,1)), log10(Ncpoints(ii,2)),'Marker','o', 'Color', 'k',...
        'MarkerFaceColor', ccol(currcol, :))
end
axis equal; xlim([3.55 4.6]); ylim([3.55 4.6]);
xticks([3.6 4.1 4.6]); yticks([3.6 4.1 4.6])
xlabel('First scan');
title('log (number of cells)')




p(2,2,2).select();cla;
line([0 2], [0 2], 'Color','k','LineStyle','--')
for ii = 1:Nmice
    currcol = (datepoints(ii,2)>50) + 1;
    line(cellsnr(ii,1),cellsnr(ii,2),'Marker','o', 'Color', 'k',...
        'MarkerFaceColor', ccol(currcol, :))
end
axis equal; xlim([0.55 1.4]); ylim([0.55 1.4])
yticks([0.6 1 1.4]); xticks([0.6 1 1.4])
title('Cell SNR')
xlabel('First scan');

%%
savepath = 'S:\ElboustaniLab\#SHARE\Documents\Dimos\ALiCe';
filenamepdf = 'manuscript_figure.pdf';
p.export(fullfile(savepath,filenamepdf), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-r300');
%%
