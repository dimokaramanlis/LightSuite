function regopts = bulkAlignToAllen(sliceinfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
regopts.howtoperm = [3 1 2];
%--------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
alignedvol = loadLargeSliceVolume(sliceinfo.slicevolfin, 1);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Loading and processing Allen Atlas template... '); tic;
allenres = 10;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv               = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
tvreg            = imresize3(tv, allenres/sliceinfo.px_register);
limskeep         = [55, size(tvreg, 1)-100]; % exclude cerebellum and olfactory
tv_cloud         = extractHighSFVolumePoints(tvreg, sliceinfo.px_register, limskeep);
Npts             = tv_cloud.Count;
tv_cloud_use     = pcdownsample(tv_cloud,'random', 10000/Npts, 'PreserveStructure',true);

fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
fprintf('Extracting points of interest from data... '); tic;

scalefilter = 100/sliceinfo.px_register;
finsize     = ceil(sliceinfo.size_proc*sliceinfo.px_process/sliceinfo.px_register);
alignedvol  = single(alignedvol);
bval        = median(single(sliceinfo.backvalues(1,:)));
alignedvol  = (alignedvol - bval)/bval;
alignedvol(alignedvol<0) = 0;
alignedvol = imresize3(alignedvol, [finsize sliceinfo.Nslices]);


% we perform spatial bandpass filtering
imhigh          = spatial_bandpass(alignedvol, scalefilter, 3, 3, sliceinfo.use_gpu);
thresuse        = quantile(imhigh(randperm(numel(imhigh), 1e5)),0.99,'all')/2;
idxcp           = find(imhigh>thresuse);
[row,col,slice] = ind2sub(size(imhigh), idxcp);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------

%%
zvals = (slice * sliceinfo.slicethickness)/sliceinfo.px_register;
xvals = col ;
yvals = row;
Xinput = [gather(yvals), gather(zvals), gather(xvals)];
% Xinput = Xinput - mean(Xinput);
pcall  = pointCloud(Xinput);
pcplot = pcdownsample(pcall, 'random', 0.05, 'PreserveStructure',true);
clf;
scatter3(pcplot.Location(:,1),pcplot.Location(:,2),pcplot.Location(:,3), 2, 'filled');
hold on;
scatter3(tv_cloud_use.Location(:,1),tv_cloud_use.Location(:,2),tv_cloud_use.Location(:,3), 2, 'filled')

%%
fprintf('Coarse alignment of slice volume to Allen atlas... '); tic;
[tformrigid, pcreg] = pcregistercpd(tv_cloud_use, pcplot, "Transform","Rigid",...
    "Verbose",true,"OutlierRatio",0.00, 'MaxIterations', 150, 'Tolerance', 1e-6);

% [r,t,data3] = icp(pcplot.Location, tv_cloud_use.Location, 200, 10, 1, 1e-6);
% tformrigid = rigidtform3d(r,t);
% pcshowpair(pcplot, pointCloud(data3'))

figure;
pcshowpair(pcplot, pcreg)

fprintf('Done! Took %2.2f s\n', toc); 
tvtrans = pctransform(tv_cloud, tformrigid);
%%
volsamp  = permute(alignedvol, regopts.howtoperm);
raatlas  = imref3d(size(tvreg),  1, 1, 1);
rasample = imref3d(size(volsamp), 1, sliceinfo.slicethickness/ sliceinfo.px_register, 1);
%%
%% adjusting the global transformation
% for each slice, we find the best +- fit and re-adjust the z-value

Nerrs    = 11;
valsuse  = linspace(-3, 3, Nerrs)*sliceinfo.slicethickness/sliceinfo.px_register;
allerrs  = nan(sliceinfo.Nslices, Nerrs, 2);

for islicecurr = 1:sliceinfo.Nslices
    fprintf('Aligning slice %d to atlas... ', islicecurr); 
    slicetic = tic;
    Xcurr = Xinput(slice == islicecurr, :);
    apcurr     = Xcurr(1,2);
    % Xcurr(:,2) = Xcurr(:,2) + randn(size(Xcurr,1),1)*sliceinfo.slicethickness/sliceinfo.px_register/4;

    for ii = 1:Nerrs
        iatlascurr  = abs(tvtrans.Location(:,2) - (apcurr+valsuse(ii))) < sliceinfo.slicethickness/sliceinfo.px_register/2;
        if nnz(iatlascurr) < 100
            continue;
        end
        pcatlascurr = pointCloud(tvtrans.Location(iatlascurr, :));
        pcatlascurr = pcdownsample(pcatlascurr, 'nonuniformGridSample',60, 'PreserveStructure',true);
        pcslicecurr = pointCloud(Xcurr);
        pcslicecurr = pcdownsample(pcslicecurr, 'nonuniformGridSample', 6, 'PreserveStructure',true);

        [R,T,data2] = icp(pcslicecurr.Location(:,[1 3]), pcatlascurr.Location(:,[1 3]), 100, 10, 1, 1e-6);
        data2 = data2';
        res(1) = mean(min(pdist2(data2, pcslicecurr.Location(:,[1 3])),[],2));
        res(2) = mean(min(pdist2(data2, pcslicecurr.Location(:,[1 3])),[],1));
        % [tformrigids, pcregs, res] = pcregistercpd(pcatlascurr, pcslicecurr, "Transform","Rigid",...
        %     "Verbose",false,"OutlierRatio",0.00, 'MaxIterations', 200, 'Tolerance', 1e-6);
        % data2 = pcregs.Location(:, [1 3]);
        allerrs(islicecurr, ii, :) = res;
        % clf;
        % subplot(1,2,1)
        % plot(pcatlascurr.Location(:,3), pcatlascurr.Location(:,1), 'r.',...
        %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
        % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
        % subplot(1,2,2)
        % plot(data2(:,2), data2(:,1), 'r.',...
        %     pcslicecurr.Location(:,3), pcslicecurr.Location(:,1), 'k.')
        % axis equal; axis tight;ax = gca; ax.YDir = 'reverse';
        % title(sprintf('Error: %3.3f, %3.3f', res(1), res(2)))
        % pause;
    end

    fprintf('Done! Min error %2.2f. Took %2.2f s\n', min(allerrs(islicecurr,:)), toc(slicetic));
    % 
    % % pcshowpair(pcslicecurr, pcregs)
    % alltrans(islicecurr,:) = tformrigids.Translation;
    % 
    % 
    % Tmat      = tformrigid.A*tformrigids.A;
    % tformcurr = affinetform3d(Tmat);
    % tvinsamp = imwarp(tvreg, raatlas, tformcurr, 'OutputView', rasample);
    % tvinsamp(islicecurr,:,:)
    
end
%%
volsamp  = permute(alignedvol, regopts.howtoperm);
raatlas  = imref3d(size(tvreg),  1, 1, 1);
rasample = imref3d(size(volsamp), 1, sliceinfo.slicethickness/ sliceinfo.px_register, 1);
tvinsamp = imwarp(tvreg, raatlas, tformrigid, 'OutputView', rasample);
%%
figure;
for ii = 1:sliceinfo.Nslices
    subplot(1,2,1); 
    imagesc(squeeze(tvinsamp(ii,:,:))); 
    axis equal; axis tight;
    subplot(1,2,2); 
    imagesc(squeeze(volsamp(ii,:,:))); 
    axis equal; axis tight;
    pause; 
end

end