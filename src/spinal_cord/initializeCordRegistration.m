function regopts = initializeCordRegistration(regopts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
samppts     = regopts.smpts;
tvpts       = regopts.tvpts;
Nslices     = range(regopts.ikeeprange) + 1;
Ntargetatlas = 1e5;
pcdownatlas  = reducePoints(tvpts, Ntargetatlas);
%==========================================================================
targetcent          = median(tvpts(:, [1 2]));
%==========================================================================
mancurdp            = fullfile(regopts.lsfolder, "spinal_alignment_opt.mat");
mandpdata           = load(mancurdp);
align_out           = mandpdata.align_out;
tforms              = computeStraighteningTransforms(align_out, targetcent, 90);
regopts.slicetforms = tforms;
%==========================================================================
% apply initial transforms
fprintf('Applying straightening transforms... '); tic;
sizetv      = size(regopts.tv, [1 2]);
raout       = imref2d(sizetv);
straightvol = tranformCordImagesSlices(regopts.regvol, tforms, raout);
fprintf('Done! Took %2.2f s.\n', toc);
%==========================================================================
%%
% p = panel();
% p.pack('v', 2);
% cimlims = getImageLimits(straightvol, 0.001);
% p(1).select();
% imagesc(squeeze(max(regopts.regvol,[],1)),cimlims)
% text(Nslices*0.95, 50, '4 mm', 'Color','w', 'HorizontalAlignment','right', 'FontSize',24)
% line(Nslices*0.95 + [-200 0], 20*[1 1], 'Color','w', 'LineWidth',2)
% 
% axis image off;  ax = gca; ax.Colormap = gray;
% p(2).select();
% imagesc(squeeze(max(straightvol,[],1)), cimlims)
% axis image off; ax = gca; ax.Colormap = gray;

%%
%==========================================================================
% clean and downsample points
thresout    = 2;
centers     = [align_out.fit_x align_out.fit_y];  
ikeepori    = removeOutliers(samppts, centers(samppts(:,3),:), thresout); 
samppts2    = tranformCordPointsSlices(samppts(ikeepori, :), tforms);
pcdownsamp  = reducePoints(samppts2, Ntargetatlas);
%==========================================================================
% move atlas to sample space
oririgid   = rigidtform3d(eye(3), [0 0 Nslices/2 - size(regopts.tv,3)/2]);
pcatlasreg = oririgid.transformPointsForward(pcdownatlas);
%==========================================================================
% let's fit the initial simiarity transform
fprintf('Initial similarity transform... '); tic;
[yreg,bfit] = pcregisterBCPD(pcatlasreg, pcdownsamp, 'TransformType','Similarity',...
    'BCPDPath', regopts.bcpdpath, 'OutlierRatio', 0.01, ...
    NormalizeCommon = true, Beta = 15, Verbose = false, ConvergenceTolerance=1e-8);
fprintf('Done! Took %2.2f s.\n', toc);
%==========================================================================
%%
% figure; hold on;
% % scatter3(yreg(:,1), yreg(:,2), yreg(:,3), 1,'filled', 'Color','k')
% % scatter3(pcdownsamp(:,1), pcdownsamp(:,2), pcdownsamp(:,3), 1,'filled', 'Color','r')
% % 
% % axis equal;
% cmapcurr = cbrewer('qual','Set1',255);
% cmapcurr = [0 0 0; cmapcurr];
% imagesc(squeeze(regopts.av(65,:,:)))
% axis image off;
% ax = gca; ax.Colormap = cmapcurr;
% 
% % imagesc(squeeze(max(regopts.av,[],1)))
%%
transinit = affinetform3d(oririgid.A*bfit.A);
refsample = imref3d(size(straightvol));
refatlas  = imref3d(size(regopts.tv));

tvtemp   = medfilt3(regopts.tv);
atlasuse = imwarp(tvtemp, refatlas, transinit, 'OutputView', refsample);
%==========================================================================
% save similarity volume for inspection
avsim  = imwarp(regopts.av, refatlas, transinit, 'nearest', 'OutputView', refsample);

volmax  = single(quantile(straightvol,0.999,'all'));
volplot = uint8(255*single(straightvol)/volmax);
cf      = plotCordAnnotation(volplot, avsim);
savepngFast(cf, regopts.lsfolder, sprintf('registration_initial_similarity'), 400, 1);
close(cf);
%==========================================================================
% perform affine registration using image information
%%
[~, ~, tformpath, ~] = performElastixAffineRegistration(atlasuse,straightvol, 1, regopts.lsfolder);
newtrans = affinetform3d(parse_elastix_tform(tformpath));

transaff = affinetform3d(transinit.A*newtrans.A);
avshow   = imwarp(regopts.av, refatlas, transaff, 'nearest', 'OutputView', refsample);
cf       = plotCordAnnotation(volplot, avshow);
savepngFast(cf, regopts.lsfolder, sprintf('registration_initial_affine'), 400, 1);
close(cf);
%%
%==========================================================================
regopts.affine_atlas_to_samp = transaff;
regopts.straightvol          = straightvol;
% save registration
save(fullfile(regopts.lsfolder, 'regopts.mat'), '-struct', 'regopts')

%==========================================================================
end

function ptsout = reducePoints(pts, Ntarget)
Nlevels = max(pts(:,3));
downfac = Ntarget/Nlevels;
ptsout  = cell(Nlevels, 1);
for ii = 1:Nlevels
    icurr      = pts(:, 3) == ii;
    Ndownatlas = double(max(6, floor(nnz(icurr)/downfac)));
    ptscurr    = pcdownsample(pointCloud(pts(icurr,:)), 'nonuniformGrid', Ndownatlas);
    ptsout{ii}  = ptscurr.Location;
end
ptsout = cat(1, ptsout{:});

% Ndownatlas   = max(6, floor(size(pts,1)/Ntarget));
% ptsout       = pcdownsample(pointCloud(pts), 'nonuniformGrid', Ndownatlas);
% ptsout       = ptsout.Location;
end

function ikeep = removeOutliers(pts, centerval, thresuse)

ikeep       = true(size(pts,1), 1);
allds       = sqrt(sum((pts(:, [1 2]) - centerval).^2,2));
Nmax        = max(pts(:, 3));
dsperslice  = accumarray(pts(:, 3), allds, [Nmax 1], @robustStd);
medperslice = accumarray(pts(:, 3), allds, [Nmax 1], @median);
medperslice = movmedian(medperslice, 50);
dsperslice  = movmedian(dsperslice,  50);
inoise      = find((allds-medperslice(pts(:,3)))./dsperslice(pts(:,3)) > thresuse);


% [~,isort]         = sort(pts(:, 3));
% inoise            = find(allds(isort)./movmedian(allds(isort), Nsmooth) > 2);
% inoise            = isort(inoise);
if numel(inoise) > 0
    fprintf('Found %d outliers around the cord\n', numel(inoise));
    ikeep             = ~ismember(1:size(pts,1), inoise)';
end
end