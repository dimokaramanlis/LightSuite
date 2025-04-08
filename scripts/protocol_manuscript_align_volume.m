mouseid      = 'YX012';
pxsize       = [1.44 1.44 5];
newdataloc   = 'F:\test_manuscript';
cdirfind     = dir(fullfile(newdataloc, sprintf('*%s*', mouseid), '**', '*.tif'));
findsavepath = dir(fullfile(newdataloc, sprintf('*%s*', mouseid)));
savepath     = fullfile(findsavepath.folder, findsavepath.name);

% percent of where match is expected

% areaslook     = [0.2 0.45; 0.3 0.7; 0.18 0.6]; % yx012 DONE
% areaslook     = [0.2 0.45; 0.3 0.7; 0.25 0.7]; % yx004 DONE
% areaslook     = [0.27 0.48; 0.28 0.72; 0.22 0.63]; % yx003 DONE
% areaslook     = [0.18 0.42; 0.33 0.67; 0.25 0.75]; % dk028 DONE
% areaslook       =[0.22 0.46; 0.31 0.69; 0.06 0.65]; % dk027 DONE
% areaslook     = [0.1 0.4; 0.25 0.75; 0.25 0.75]; % dk025 DONE
% areaslook     = [0.2 0.52; 0.32 0.68; 0.05 0.55]; % dk026 DONE
% areaslook     = [0.23 0.45; 0.3 0.7; 0.2 0.6]; % yx011 DONE

%%
%=========================================================================
% FIRST ANALYSIS PART
%=========================================================================
res           = struct();
res.pxsize    = pxsize;
% areaslook = [0.1 0.45; 0.3 0.9; 0.1 0.75]; % dk025 looks good, rotation?
%=========================================================================
% new imaging
newdate      = datetime(cdirfind(1).date);
res.newdate  = newdate;
%=========================================================================
% oldimaging
dpfind          = fullfile('F:\imaging', sprintf('%s*', mouseid), '**','*.tif');
allfilesfind    = dir(dpfind);
[allun, ~, iun] = unique({allfilesfind(:).folder}');
[~, ifolder]    = max(accumarray(iun, 1, [numel(allun) 1], @sum));
origdatafolder  = allun{ifolder};
dpolddate       = fullfile('F:\imaging', sprintf('%s*', mouseid), '**','*.ini');
dpoldfile       = dir(dpolddate);
res.olddate     = datetime(dpoldfile.date);
%=========================================================================
% extract the new volume

newfilelist = fullfile({cdirfind(:).folder}', {cdirfind(:).name}');
Nuse        = 500;
Nfileskeep  = min(numel(newfilelist), Nuse);
tfile       = imread(newfilelist{1});
[Ny, Nx]    = size(tfile);
ikeep       = round(numel(newfilelist)/2 + (-Nfileskeep/2:Nfileskeep/2));
ikeep(ikeep<1 | ikeep>numel(newfilelist)) = [];
Nz          = numel(ikeep);

voltot    = zeros(Ny, Nx, Nz, 'uint16');
signalvol = zeros(Ny, Nx, Nz, 'uint16');

diaminpx     = 14./1.44;
opts.maxdff  = 12;
normfac      = single(intmax("uint16"))/opts.maxdff;
medwithfull  = 2*ceil((2.5*diaminpx)/2) + 1;

tic; msg = [];
for ifile = 1:Nz
    currim = imread(newfilelist{ikeep(ifile)});
    voltot(:, :, ifile) = medfilt2(currim, [3 3]);

    backframe = single(medfilt2(currim, [medwithfull medwithfull]));
    cellimage = (single(voltot(:, :, ifile)) -backframe)./backframe;
    cellimage(isinf(cellimage) | isnan(cellimage)) = 0;
    signalvol(:,:,ifile) = uint16(cellimage*normfac);


    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Files %d/%d. Time elapsed %2.2f s...\n', ifile, Nz,toc);
    fprintf(msg);
end

rng(1);
isample              = randperm(numel(voltot), 2e6);
bottomval            = quantile(voltot(isample), 0.001, 'all');
topval               = quantile(voltot(isample), 0.999, 'all');
DR                   = topval-bottomval;
res.newhistogramvals = [bottomval topval];
res.filerangenewvol  = ikeep;
fprintf('Dynamic range of %s new is %d\n', mouseid, DR)
%=========================================================================
% extract cells from the new volume
[opts.Ny, opts.Nx, opts.Nz] = size(signalvol);
newcelllocs = extractCellsFromVolume(signalvol, opts);
% newcelllocs(:,3) = interp1(1:numel(ikeep), ikeep, newcelllocs(:,3));
% let's save the cell locations
res.newcelllocs = newcelllocs;

samplepath = fullfile(savepath, sprintf('signalvolnew_%s.tif', mouseid));
options.compress = 'lzw';
options.message  = false;
if exist(samplepath, 'file')
    delete(samplepath);
end
saveastiff(voltot, samplepath, options);

res.newsignalvol = samplepath;
save(fullfile(savepath, 'resnew.mat'),'-struct', 'res')
%%
%=========================================================================
% SECOND ANALYSIS PART
%=========================================================================

% lets's extract cells from the original volume in a fast way
[voloriginsig, voloriginact, yxzranges] = preprocessOriginalVolume(origdatafolder, areaslook);
[opts.Ny, opts.Nx, opts.Nz] = size(voloriginsig);
opts.maxdff = 12;
cell_locs_old = extractCellsFromVolume(voloriginsig, opts);
res.yxzrangesold = yxzranges;
clear voloriginsig;
%=========================================================================
%%
% algorithm for matching points
iusenew   = newcelllocs(:,4)>8;
iuseold   = cell_locs_old(:,4)>8;

oldpoints = cell_locs_old(iuseold, 1:3).* pxsize;
newpoints = newcelllocs(iusenew, 1:3)  .* pxsize;

rng(1);
origpc     = pointCloud(oldpoints);
% origpcdown = origpc;
origpcdown = pcdownsample(origpc, "nonuniformGridSample",15, 'PreserveStructure', true);
newpc      = pointCloud(newpoints);
% newpcdown  = newpc;
newpcdown  = pcdownsample(newpc, "nonuniformGridSample", 15, 'PreserveStructure', true);


[tform, movingRegdown] =  pcregistercpd(origpcdown,newpcdown,'Transform', 'Rigid',...
    'Verbose',true,'MaxIterations',200,'Tolerance',1e-7, 'OutlierRatio', 0.0);
movingReg = pctransform(origpc, tform);

figure;
subplot(1,2,1)
pcshowpair(origpc,newpc,'MarkerSize',25)
subplot(1,2,2)
pcshowpair(movingReg,newpc,'MarkerSize',25)
tform_combined = tform;



%%
voltnewref = imref3d(size(voltot), 1.44, 1.44, 5);
voltoldref = imref3d(size(voloriginact), 1.44, 1.44, 5);

volfin = imwarp(voloriginact, voltoldref, tform_combined, 'linear','OutputView', voltnewref);

isample   = randperm(numel(volfin), 2e6);
bottomvalold = quantile(volfin(isample), 0.001, 'all');
topvalold    = quantile(volfin(isample), 0.999, 'all');
DRold =  topvalold-bottomvalold;
fprintf('Dynamic range of %s old is %d\n', mouseid, DRold)

%%
%390-392
rangeshow = 100:120; %380-400 for DK026
figure; p = panel(); p.pack('h',2); 
p.de.margin = 1;
p(1).select();cla;
imagesc(max(volfin(:,:,rangeshow),[],3),[150 4000]); 
ax = gca; ax.Colormap = copper; axis equal; axis tight;
ax.Visible ='off'; ax.Title.Visible = 'on';
title(sprintf('%s before', mouseid))
p(2).select();cla;
imagesc(max(voltot(:,:,rangeshow),[],3),[50 4000])
ax = gca; ax.Colormap = copper; axis equal; axis tight;
ax.Visible ='off'; ax.Title.Visible = 'on';
title(sprintf('%s after', mouseid))
%%
oldcellsmapped = tform_combined.transformPointsForward(cell_locs_old(:,1:3).* pxsize);
ikeepoldx = oldcellsmapped(:,1)>voltnewref.XWorldLimits(1) & oldcellsmapped(:,1)<voltnewref.XWorldLimits(2);
ikeepoldy = oldcellsmapped(:,2)>voltnewref.YWorldLimits(1) & oldcellsmapped(:,2)<voltnewref.YWorldLimits(2);
ikeepoldz = oldcellsmapped(:,3)>voltnewref.ZWorldLimits(1) & oldcellsmapped(:,3)<voltnewref.ZWorldLimits(2);
ikeepold  = ikeepoldx & ikeepoldz& ikeepoldy;
oldcellsmapped = oldcellsmapped(ikeepold, :);

res.areaslook = areaslook;
res.oldhistogramvals = [bottomvalold topvalold];
res.oldcelllocations = [oldcellsmapped cell_locs_old(ikeepold, 4:5)];

samplepathold = fullfile(savepath, sprintf('signalvolold_%s.tif', mouseid));
if exist(samplepathold, 'file')
    delete(samplepathold);
end
saveastiff(volfin, samplepathold, options);

res.oldsignalvol = samplepathold;
save(fullfile(savepath, 'resfull.mat'),'-struct', 'res')

