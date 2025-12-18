function regopts = prepareCordAtlasForRegistration(cordvol, opts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
fprintf('Loading cord atlas and extracting point clouds... '); tic;
[tv, av, tvpts]     = loadSpinalCordAtlas(opts.dpspineatlas, opts.sampleres);
fprintf('Done! Took %2.2f s.\n', toc)
%--------------------------------------------------------------------------
iregchan  = opts.regchan;
if opts.Nchan == 1
    iregchan = 1;
end
assert(iregchan<=opts.Nchan , 'Your registration channel is out of range');
%--------------------------------------------------------------------------
% we assume the long dimension is first
regvol           = cordvol(:, :, :, iregchan);
[~, im]          = max(opts.orisize);
axperm           = setdiff(1:3, im);
opts.sampleperm  = [axperm im];
regvol           = permute(regvol, opts.sampleperm);
regvolsize       = size(regvol);
fprintf('Found the long axis in position %d, moved to last\n', im)
%--------------------------------------------------------------------------
fprintf('Segmenting the cord... '); tic;
seuse             = strel('cuboid', [3 3 3]);
isamprand         = randperm(numel(regvol), 2e4);
sampsuse          = regvol(isamprand);
backval           = mode(sampsuse(sampsuse>0));
regvol(regvol==0) = backval;
binvol            = imbinarize(regvol);
binvol            = imdilate(binvol, seuse);
cc                = bwconncomp(binvol, 18);
cinfo             = regionprops3(cc, 'Volume', 'VoxelList');
newvol            = false(size(binvol));
[~, imax]         = max(cinfo.Volume);
voxset            = cinfo.VoxelList{imax};
idsshow           = sub2ind(size(newvol), voxset(:,2), voxset(:,1), voxset(:,3));
newvol(idsshow)   = true;
fprintf('Done! Took %2.2f s.\n', toc);
%--------------------------------------------------------------------------
% we then find whether rostrocaudal works
cordarea          = squeeze(sum(newvol, 1:2));
Ncheck            = ceil(0.05*numel(cordarea));
frontsum          = mean(cordarea(1:Ncheck));
backsum           = mean(cordarea(end-Ncheck+1:end));
opts.tofliprc     = frontsum < backsum;
if opts.tofliprc
    fprintf('Spinal cord is in caudorostral direction, flipping...\n')
    tv   = flip(tv, 3);
    av   = flip(tv, 3);
    tvpts(:, 3) = size(tv, 3) - tvpts(:, 3);
end
opts.cordarea = cordarea;
%--------------------------------------------------------------------------
% find and remove brain parts
ihigh    = cordarea > (median(cordarea) + 3*robustStd(cordarea));
fprintf('%2.2f%% of the cord has larger size than expected, skipping\n', mean(ihigh)*100)
ikeep    = [1 size(regvol, 3)];
if mean(ihigh) > 0.01
    if opts.tofliprc 
        ilast = find(ihigh==1, 1);
        ikeep = [1 ilast+1];
    else
        ifirst = find(ihigh==0, 1);
        ikeep  = [ifirst-1 size(regvol, 3)];
    end
end
opts.ikeeprange = ikeep;
%--------------------------------------------------------------------------
% Example Data (Replace with your actual 'area' vectors)
% Sample: The fixed "target" (e.g., your spinal cord sample)
% Atlas: The "moving" reference you want to stretch/warp
binatlas  = imbinarize(tv);
atlasarea = squeeze(sum(binatlas, 1:2));
cordmatch = cordarea(ikeep(1):ikeep(2));
%%
sample = cordmatch./median(cordmatch(cordmatch>0)); 
atlas  = atlasarea./median(atlasarea);     

atlas  = [0; atlas];

% 1. Run Standard DTW
% ix corresponds to sample indices, iy to atlas indices
[dist, ix, iy] = dtw(sample, atlas, 100);

% 2. Calculate the Warping Function
% We need to know: for every point in 'sample', what is the corresponding 
% index in 'atlas'?
% Since 'ix' (sample indices) will have duplicates (repeats), we average 
% the corresponding 'iy' values to find the center of the match.

len_sample = length(sample);
warp_path = zeros(1, len_sample);

for k = 1:len_sample
    % Find all indices in the warping path that correspond to sample index 'k'
    matches = (ix == k);
    
    % Take the average atlas index that matches this sample point
    warp_path(k) = mean(iy(matches)); 
end
opts.atlaswarppath = warp_path;
% % 3. Warp the Atlas to the Sample's Grid
% % We now interpolate the Atlas values at these new mapped coordinates.
% % 'linear' or 'pchip' (shape-preserving) are good options.
% atlas_aligned = interp1(1:length(atlas), atlas, warp_path, 'pchip');
% 
% % --- Visualization ---
% figure;
% subplot(3,1,1); plot(sample, 'b'); title('Fixed Sample'); grid on;
% subplot(3,1,2); plot(atlas, 'r'); title('Original Atlas'); grid on;
% subplot(3,1,3); 
% plot(sample, 'b', 'LineWidth', 2); hold on;
% plot(atlas_aligned, 'r--', 'LineWidth', 2);
% title('Atlas Warped to Match Sample'); 
% legend('Sample', 'Aligned Atlas'); grid on;
%%
%--------------------------------------------------------------------------
% reducing volume

xrange   = [max(min(voxset(:, 1)) - 1, 1) min(max(voxset(:, 1)) + 1, regvolsize(2))];
yrange   = [max(min(voxset(:, 2)) - 1, 1) min(max(voxset(:, 2)) + 1, regvolsize(1))];
revoluse  = regvol;
revoluse  = revoluse(yrange(1):yrange(2), xrange(1):xrange(2), :);
newvol    = newvol(yrange(1):yrange(2), xrange(1):xrange(2), :);
%--------------------------------------------------------------------------
fprintf('Extracting point clouds from sample... '); tic;
ptvol        = spinalCordPointCloud(revoluse, ikeep, newvol);
fprintf('Done! Took %2.2f s.\n', toc);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
opts.tv     = tv;
opts.av     = av;
opts.xrange = xrange;
opts.yrange = yrange;
opts.regvol = revoluse(:, :, ikeep(1):ikeep(2));
opts.tvpts  = tvpts;
opts.smpts  = ptvol;
regopts     = opts;
dpsave      = fullfile(opts.lsfolder, 'regopts.mat');
save(dpsave, '-struct', 'regopts')
%--------------------------------------------------------------------------
end



% binvol = false(size(regvol));
% Nslices   = opts.orisize(1);
% Nbatch    = 4;
% batchsize = ceil(Nslices/Nbatch);
% for ibatch = 1:Nbatch
%     istart = (ibatch-1)*batchsize + 1;
%     iend   = min(ibatch*batchsize, Nslices);
%     binvol(istart:iend, :, :) = imbinarize(regvol(istart:iend, :, :));
% end