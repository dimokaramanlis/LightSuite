function [flipdecisions, tformarray] = alignConsecutiveSlicesFlip(pointcell, indsmatch, imwidth)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nslices       = numel(indsmatch);
flipdecisions = false(Nslices, 1);
% tformarray    = rigidtform2d; % initialize as identity
tformarray    = affinetform2d; % initialize as identity

for islice = 1:numel(indsmatch)-1
    ifixed    = indsmatch(islice);
    imoving   = indsmatch(islice + 1);
    fprintf('Aligning slice %d (moving) to %d (fixed)... ', imoving, ifixed); 
    slicetic = tic;
    movingcloud = pointcell(imoving, :);
    fixedcloud  = pointcell{ifixed, 1};

    errscurr  = nan(2,1);
    pcregcurr = cell(2,1);
    tformcurr = cell(2,1);
    for iside = 1:2
        [tformcurr{iside, 1}, pcregcurr{iside, 1}, errscurr(iside, 1)] = pcregistercpd(...
            movingcloud{iside}, fixedcloud, ...
            "Transform","Rigid",'Verbose',false,...
            'Tolerance',1e-7,'OutlierRatio',0.00,'MaxIterations', 100);
    end
    % for iside = 1:2
    %     [tformcurr{iside, 2}, pcregcurr{iside, 2}, errscurr(iside, 2)] = pcregistercpd(...
    %         fixedcloud, movingcloud{iside}, ...
    %         "Transform","Rigid",'Verbose',false,...
    %         'Tolerance',1e-6,'OutlierRatio',0.00,'MaxIterations', 100);
    % end
    % errscurr
    % [~, imin]= min(errscurr,[],'all');
    % errscurr   = harmmean(errscurr, 2);
    % errsdecide   = min(errscurr, [], 2);
    errorratio  = (errscurr(2) - errscurr(1))./errscurr(1);
    % [toflip, ~] = ind2sub([2 2], imin);

    if errorratio < 0
        stradd = ' Flipping!';
        flipdecisions(islice + 1) = true;
        tformuse = rigid3dToRigid2d(tformcurr{2});
        M_flip_x = [-1  0  imwidth; 0 1  0; 0  0  1];
    else
        tformuse = rigid3dToRigid2d(tformcurr{1});
        stradd = '';
        M_flip_x = [1  0  0; 0 1  0; 0  0  1];
    end
    % tformarray(islice+1) = tformuse;
    M_new                = tformuse.A * M_flip_x;
    tformarray(islice+1) = affinetform2d(M_new);

    fprintf('Done! Err change %2.2f%%.%s Took %2.2f s\n', errorratio*100, stradd, toc(slicetic));
    %--------------------------------------------------------------------------
    % % for debugging
    clf;
    subplot(1,3,1)
    plot(movingcloud{1}.Location(:,1), movingcloud{1}.Location(:,2), 'r.',...
        fixedcloud.Location(:,1), fixedcloud.Location(:,2), 'k.')
    axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    subplot(1,3,2)
    plot(pcregcurr{1}.Location(:,1), pcregcurr{1}.Location(:,2), 'r.',...
        fixedcloud.Location(:,1), fixedcloud.Location(:,2), 'k.')
    axis equal; axis tight;ax = gca; ax.YDir = 'reverse';
    title(sprintf('Error: %3.3f', errscurr(1)))

    subplot(1,3,3)
     plot(pcregcurr{2}.Location(:,1), pcregcurr{2}.Location(:,2), 'r.',...
        fixedcloud.Location(:,1), fixedcloud.Location(:,2), 'k.')
    axis equal; axis tight;ax = gca; ax.YDir = 'reverse';
    title(sprintf('Flipped error: %3.3f', errscurr(2)))
    %--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
% figure;
% subplot(1,3,1)
% pcshowpair(allclouds{itemplate-1,1}, allclouds{itemplate, 1});
% view(2); ax = gca; ax.YDir = 'reverse';
% subplot(1,3,2)
% pcshowpair(pcregcurr{1}, allclouds{itemplate, 1})
% view(2); ax = gca; ax.YDir = 'reverse';
% title(sprintf('Error: %3.3f', errscurr(1)))
% 
% subplot(1,3,3)
% pcshowpair(pcregcurr{2}, allclouds{itemplate, 1})
% view(2); ax = gca; ax.YDir = 'reverse';
% title(sprintf('Flipped error: %3.3f', errscurr(2)))
% % fprintf('Normal err = %3.3f, flipped err = %3.3f\n', err, errflip)
% 
% if errscurr(1)<errscurr(2)
%     tformuse = rigid3dToRigid2d(tformcurr{1});
%     imuse    = volregister(:, :, itemplate-1);
% else
%     tformuse = rigid3dToRigid2d(tformcurr{2});
%     imuse    = flip(volregister(:, :, itemplate-1), 2);
% end
% 
% raref        = imref2d(size(imuse));
% imregistered = imwarp(imuse, tformuse, 'linear','OutputView',raref);
% figure;
% imshowpair(volregister(:, :, itemplate), imregistered)
%--------------------------------------------------------------------------
end