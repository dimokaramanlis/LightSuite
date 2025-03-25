%--------------------------------------------------------------------------
% this is the stitched file path

dp = 'S:\ElboustaniLab\#SHARE\Data\DK026\Anatomy\RES(641x534x87)\000000\000000_0-1450';
tiffpaths   = dir(fullfile(dp, '*.tif'));

slicesplot = [70 40 25];
facdown = 16;
xyres   = 1.44;
dispxy  = [15000 000];

clf; hold on;
for ii = 1:numel(slicesplot)
    islice = slicesplot(ii);
    currim = imread(fullfile(tiffpaths(islice).folder, tiffpaths(islice).name));
    xx = (1:size(currim, 2))*xyres * facdown;
    xx = xx + (ii-1)*dispxy(1);
    yy = (1:size(currim, 1))*xyres * facdown;
    yy = yy + (ii-1)*dispxy(2);
    imagesc(xx, yy, currim, [0 2500])
    if ii == 3
        text(min(xx) + 200, min(yy)+1000, '2 mm', 'FontSize',22, 'Color','w')
        line(min(xx) + 200 + [0 2000], min(yy) + [1 1]* 2000, 'Color', 'w',...
            'LineWidth', 2)
    end
end
ax = gca; ax.YDir = "reverse"; ax.Colormap = gray;
axis equal; axis tight; ax.Visible = 'off';