% set main data path - where the tiffs are
dpspinesample   = 'D:\spine_registration\sample2';
%% load sample and atlas - set resolution and channel for registration
sampleres          = 20; % in micrometers
[cordvol, opts]    = readSpinalCordSample(dpspinesample, sampleres);
opts.dpspineatlas  = 'D:\AllenAtlas\extra_spine\allen_cord_20um_v1.1';
opts.regchan       = 2; % choose registration channel
regopts            = prepareCordAtlasForRegistration(cordvol, opts);
%% (manual) trace spinal cord to untwist and unbend
regopts = loadRegOpts(dpspinesample);
spinal_cord_aligner(regopts);

%% (auto) initialize registration
sctransformparams = initializeCordRegistration(regopts);



%%
tformslices(ilast, 1) = rigidtform2d;

globalSpinalCordFit()

transformtypes = {'rigid', 'affine'};
transformsteps = [1 1 1];
errall = nan(numel(transformsteps), 1);
for istep = 1:numel(transformsteps)
    currtranstype               = transformtypes{transformsteps(istep)};
    fprintf('Optimization step %d/%d: %s\n', istep, numel(transformsteps), currtranstype)
    [tformrigid, errall(istep)] = alignCordAtlasToSample(tv_cloud_use, ls_cloud_use, tformslices);
    tformslices                 = refineSampleFromAtlas(tv_cloud, pcsample, tformrigid, currtranstype);
    [rrx, rry, rrz]             = reportRotationAngles(tformrigid.R);
    fprintf('%s\n', repmat('=', [1 75]));
end

%%
[Nslices, Ny, Nx] = size(volcurr);
batchsize  = 100;
Nbatches   = ceil(Nslices/batchsize);
[XX, YY]  = meshgrid(1:Nx, 1:Ny);
xybuffer  = 50;
%%
fprintf('Straightening the cord... '); tic;
centsxz   = nan(Nslices, 2);
boxranges = nan(Nbatches, 2);
for ibatch = 1:Nbatches
    % we fit a tube model per batch
    istart   = (ibatch - 1) * batchsize + 1;
    iend     = min(ibatch * batchsize, Nslices);

    iptscurr   = ls_cloud.Location(:, 2) > istart & ls_cloud.Location(:, 2) < iend;
    ptcurr     = pointCloud(ls_cloud.Location(iptscurr, :));
    volbatch   = volcurr(istart:iend, :, :);
    medfit     = squeeze(quantile(volbatch,0.99,1));
    medbin     = imbinarize(medfit);
    [row, col] = find(medbin);
    miny       = max(quantile(row, 0.01) - xybuffer, 1);
    maxy       = min(quantile(row, 0.99) + xybuffer, Ny);
    minx       = max(quantile(col, 0.01) - xybuffer, 1);
    maxx       = min(quantile(col, 0.99) + xybuffer, Nx);

    volbin = imbinarize(volbatch);
    idxfore = find(volbin);
    [yy, xx, zz] = ind2sub(size(volbatch), idxfore);
    xmed = accumarray(yy, xx, [size(volbin,1) 1], @median);
    zmed = accumarray(yy, zz, [size(volbin,1) 1], @median);
    centsxz(istart:iend, :) = [xmed zmed];
    boxranges(ibatch, :) = [ceil(maxx-minx) ceil(maxy-miny)];
    % [params] = fitgaussrf(1:Nx, 1:Ny, double(medfit));
    % % cel = getEllipseFromNewParams(params, 2);
    % clf;
    % imagesc(medfit)
    % rectangle('Position',[minx miny maxx-minx maxy - miny])
    % % line(cel(1,:), cel(2, :), 'Color', 'r')
    % pause;
    % 
    % spx      = minx:maxx;
    % spy      = miny:maxy;
    % [xx, yy] = meshgrid(spx, spy);
    % 
    % pout = fitSpineCentroidModel(spx, spy, volbin(:, spy, spx)); 
    % zmed = pout(1:size(volbin,1));
    % xmed = pout(size(volbin,1)+1:end);
    % pout = fitSpineEllipseModel(spx, spy, volbin(:, spy, spx));


    clf;hold on;
    scatter3(ptcurr.Location(:,1),ptcurr.Location(:,2),ptcurr.Location(:,3),1,'filled','k')
    scatter3(xmed, istart:iend, zmed, 2, 'r')
    pause;
end
