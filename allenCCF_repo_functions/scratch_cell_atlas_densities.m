dp = 'D:\AllenAtlas\2020_annotation\ero_et_al_2018_densities.CSV';
densfile = readtable(dp);

ilayer6a  = find(contains(lower(densfile.Regions),'layer 6a'));
ilayer6b  = find(contains(lower(densfile.Regions),'layer 6b'));
ilayer5   = find(contains(lower(densfile.Regions),'layer 5'));
ilayer23  = find(contains(lower(densfile.Regions),'layer 2/3'));
ilayer4   = find(contains(lower(densfile.Regions),'layer 4'));
ilayer1   = find(contains(lower(densfile.Regions),'layer 1'));
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
alldensities = densfile.Neurons_mm_3_;
% 
% normfacs     = sum(fullsignal(idslayers,:));
% alldensities = fullsignal./normfacs;

maxc         = quantile(alldensities(cat(1, layercell{:}), :), 0.999,'all');
minc         = quantile(alldensities(cat(1, layercell{:}), :), 0.001,'all');

for ii = 1:numel(layercell)
    il = layercell{ii}; 
    p(ii).select();
    [~, isort] = sort(densfile.Regions(il));
    currdensities = alldensities(il(isort), :);
    imagesc(currdensities,[minc maxc])
    yticks(1:numel(il))
    yticklabels(densfile.Regions(il(isort)))
    ax = gca; ax.YDir = 'reverse'; ax.Colormap = magma;
    axis tight;
    title(layertitles{ii})
    xlabel('Mice'); 
end