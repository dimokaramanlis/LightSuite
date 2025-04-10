% 

datafolder = 'S:\ElboustaniLab\#SHARE\Processed\multisensory_integration\anatomy_lightsheet\data\colm';
xmlpath    = dir(fullfile(datafolder, '*561.xml'));
tileinfo   = loadColmTileInfo(fullfile(xmlpath.folder, xmlpath.name));
Ntiles     = numel(tileinfo);
procpath   = 'L:\DATA_sorted';
makeNewDir(procpath)

for itile = 1:Ntiles
    tilefiles  = dir(fullfile(datafolder, 'VW0', sprintf('LOC%03d', itile-1), '*CHN01*.tif'));
    tilefiles  = fullfile(tilefiles(1).folder, {tilefiles(:).name}');
    
    pathload   = fullfile(procpath, sprintf('LOC%03d', itile-1));
    runJetrawDecompress(pathload, tilefiles);
end





%%
% start with preprocessing
dpim     = 'L:\colm_examples';
savepath = 'C:\DATA_sorted';
[backvol, folderproc] = preprocessColmVolume(dpim, savepath);
fbin = fullfile(savepath, 'fullbrain.dat');
expdata.fbinary = fbin;
% expdata.Ny      = Ny;
% expdata.Nx      = Nx;
% expdata.Nslices = Nslices;
expdata.Ny      = 11672;
expdata.Nx      = 8464;
expdata.Nslices = 1581;
peakvalsextract = preprocessColmVolumeBatches(expdata);


%%

% now try to segment

procpath = 'C:\DATA_sorted\colm_examples';
tim     = dir(fullfile(procpath, '*.tif'));
%%

Nfiles = numel(tim);
tfile  = imread(fullfile(tim(1).folder, tim(1).name));
[Ny, Nx] = size(tfile);
tproc    = gpuArray.zeros(Ny, Nx, 'single');

peakvals = zeros(1e6, 4);
i0 = 0;
msg = []; tic;

for ifile = 1:Nfiles

    % collect all possible local maxima
    tproc(:)  = imread(fullfile(tim(ifile).folder, tim(ifile).name));
    tproc  = tproc*8/255;
    
    smin  = my_min(-tproc, [10 10], [1 2]);
    
    peaks = single(-tproc<smin+1e-3 & tproc>1.5);
    [row, col, mu] = find(peaks.*tproc);
    
    if i0+numel(row)>size(peakvals,1)
        peakvals(1e6 + size(peakvals,1), 1) = 0;
    end
    peakvals(i0 + (1:numel(row)), :) = [row col ifile*ones(size(row)) mu];
    i0 = i0 + numel(row);

    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d. Ncells %d. Time elapsed %2.2f s...\n', ifile, Nfiles, i0, toc);
    fprintf(msg);
end

peakvals = peakvals(1:i0, :);

% imagesc(tproc,[0.1 1]); axis equal; axis tight;
% line(col,row,'Marker','o','Color','r','LineStyle','none')
% axis equal; xlim([600 1100]); ylim([5600 6500])
%%
% extract 3d maxima over batches to generate candidate list
% running over batches again
% for each batch (+ overlap a bit), get brightest point, fit 3D Gaussian,
% remove point from list and all points that are within 2sigma of that
% Gaussian, and remove Gaussian from volume




%%
tim     = dir(fullfile(dpim, '*.tif'));
Nslices  = numel(tim);
tfile    = imread(fullfile(tim(1).folder, tim(1).name));
[Ny, Nx] = size(tfile);
Nside    = 2;
dslice   = -Nside:Nside;
wtsuse   = gausslinefun([0 1 1], dslice);
wtsuse   = wtsuse/sum(wtsuse);
%%
voluse   = nan(Ny, Nx, 2*Nside+1, 'single');
meduse   = nan(Ny, Nx, 2*Nside+1, 'single');

Nfilt    = [21 21];
msg = []; tic;
for islice = 1:Nslices+Nside
    voluse = circshift(voluse, -1, 3);
    meduse = circshift(meduse, -1, 3);
    meduse(:, :, end) = nan;
    if islice <= Nslices
        currim = imread(fullfile(tim(islice).folder, tim(islice).name));
        filtim = medfilt2(currim, Nfilt);
        voluse(:, :, end) = currim;
        meduse(:, :, end) = currim;
    end

    if islice >= Nside + 1
        indscrr = islice + (-Nside:Nside);
        iuse    = indscrr > 0 & indscrr<=Nslices;

        wtscurr = wtsuse(iuse);
        wtscurr = wtscurr/sum(wtscurr);
        backsignal = reshape(meduse(:,:,iuse), Ny*Nx, nnz(iuse)) * wtscurr;
        backsignal = reshape(backsignal, [Ny, Nx]);

        cellimage  = (voluse(:, :, Nside+1) - backsignal)./backsignal;
        imagesc(backsignal, [0 1e3])
        title(islice - Nside)
        drawnow
    end

    
    if mod(islice, 2) == 1 || islice == Nslices
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', islice, Nslices,toc);
        fprintf(msg);
    end
end



%%
tic;
ii = 154;
Nside = 2;
dtload = -Nside:Nside;
currim = zeros(11672, 8464, 2*Nside+1, 'uint16');
for il = 1: 2*Nside+1
    iload = dtload(il) + ii;
    currim(:,:,il) = imread(fullfile(tim(iload).folder, tim(iload).name));
end
toc;
%%
pxsize   = [1.4 1.4 5];
celldiam = 14;
% Nfilt    = 2*floor(2*celldiam./pxsize/2) + 1;
Nfilt    = [21 21 3];

tic;
backsignal = medfilt3(currim, Nfilt);
toc;
backsignal = single(backsignal);
newimage   = (single(currim) - backsignal)./backsignal;
toc;

tic;
singim = single(currim);
backsignal = medfilt2(singim, [21 21]);
newimage   = (singim - backsignal)./backsignal;
toc;

% backsignal = imgaussfilt(currim, 21);

% imagesc(, [0 5])

%%
% for background imadjust(imflatfield(backsignal, 1000))

%%
bfim = BioformatsImage(dpim);
Nz   = bfim.sizeZ;
Nx   = bfim.width;
Ny   = bfim.height;
Vdata = zeros(Ny, Nx, Nz, 'uint16');
msg = []; tic;
for ii = 1:Nz
    Vdata(:, :, ii) = bfim.getPlane(ii,1,1,1);
    if mod(ii, 20)==1 || ii ==Nz
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', ii, Nz,toc);
        fprintf(msg);
    end
end

%%
% find local maxima that are at maximum (1/pxsize)^ dimension 
% sort by intensity in the center 10 px
% do a pca and start search along the main component
% after a few samples, start rotating the component
celldiam = 12;
dpx = [8.2 8.2 5];
nfac = 1;
qthres = 0.05;
irand  = randperm(numel(Vdata),10000);
thres    = quantile(Vdata(irand),qthres,'all','method','appr');
threstop = quantile(Vdata(irand),0.9,'all');
%%
% quick thres caclulation
Nnoise = 40;
datagpu = gpuArray.zeros(Ny, Nx, 'single');
inoise = round(linspace(Npts*0.2,Npts*0.8, Nnoise));
allpts  = nan(Ny*Nx, numel(inoise), 'single');
Nfilt = 4*ceil(celldiam./dpx(1:2)) + 1;
madfac   = 1.4826;
for ii = 1:numel(inoise)
   
    idx = inoise(ii);
    datagpu(:) = volumeIdtoImage(Vdata, [idx idimuse]);
    datagpu(datagpu<thres) = 0;
    background = medfilt2(datagpu, Nfilt);
    datagpu   = datagpu - background;
    localmad  = medfilt2(abs(datagpu), Nfilt);
    datagpu   = datagpu./(localmad*madfac);
    allpts(:, ii)  = reshape(datagpu, Ny*Nx, 1);
end
threstop = quantile(allpts, 0.99,'all') * 2;

%%
idimuse = 3;
Npts = size(Vdata, idimuse);
datagpu = gpuArray.zeros(Ny, Nx, 'single');

Nwin = 4;
[xx,yy] = meshgrid(-Nwin:Nwin,-Nwin:Nwin);
tic; msg = [];
for ii = 1:Npts
   
    datagpu(:) = volumeIdtoImage(Vdata, [ii idimuse]);
    datagpu(datagpu<thres) = 0;
    imfin = bandpass_cell_filter(datagpu, 3);
    % 
    % background = medfilt2(datagpu, Nfilt);
    % datagpu   = datagpu - background;
    % localmad  = medfilt2(abs(datagpu), Nfilt);
    % datagpu   = datagpu./(localmad*madfac);
    % datagpu(localmad<1e-3) = 0;
    % imagesc(datagpu - background)
   
    
    % only keep filters centered at their max

    smin = my_min(-imfin, [2 2], [1 2]);
    peaks = single(-imfin<smin+1e-3 & -imfin<-200);
    [row, col, mu] = find(peaks);
    irem = (row<Nwin+1 | row>Ny-Nwin-1) | (col<Nwin+1 | col>Nx-Nwin-1);

    indfilts  = sub2ind([Ny Nx], yy(:) + row(~irem)', xx(:) + col(~irem)');
    allfilters = datagpu(indfilts);

    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', ii, Npts,toc);
    fprintf(msg);
end
dimiter = 3;
Ntiles = 4;
overlap = 0.1;
tiledgesx = round(linspace(1, Nx, Ntiles + 1));
tilelenx = overlap*Nx/Ntiles;
tedgesx = round(tiledgesx' + tilelenx*[-1 1]); 
tedgesx = [tedgesx(1:Ntiles,1) tedgesx(2:Ntiles+1,2)];
tedgesx(:,1) = max(tedgesx(:,1), 1);
tedgesx(:,2) = min(tedgesx(:,2), Nx);

tiledgesy = round(linspace(1, Ny, Ntiles + 1));
tileleny = overlap*Ny/Ntiles;
tedgesy = round(tiledgesy' + tileleny*[-1 1]); 
tedgesy = [tedgesy(1:Ntiles,1) tedgesy(2:Ntiles+1,2)];
tedgesy(:,1) = max(tedgesy(:,1), 1);
tedgesy(:,2) = min(tedgesy(:,2), Nx);
%%

for iy = 1:Ntiles
    indy =  tedgesy(iy, 1):tedgesy(iy, 2);
    for ix = 1:Ntiles
        indx =tedgesx(ix, 1):tedgesx(ix, 2);

        data = Vdata(indy, indx, :);
        data = single(gpuArray(data));
        % medfilt3(data, )

    end
end


%%
% go batch by batch and use the gpu...
for idim = 1:3
    
    Nmax = round(nfac*Npts*dpx(idim)/celldiam);
    
    islocalmax(Vdata, idim, 'MinSeparation', Nmin,'MaxNumExtrema', 2)
end