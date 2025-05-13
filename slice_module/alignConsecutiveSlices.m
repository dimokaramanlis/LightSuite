function tformarray = alignConsecutiveSlices(pointcell, indsmatch)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nslices       = numel(indsmatch);
tformarray    = affinetform2d; % initialize as identity

for islice = 1:numel(indsmatch)-1
    ifixed    = indsmatch(islice);
    imoving   = indsmatch(islice + 1);

    fprintf('Aligning slice %d (moving) to %d (fixed)... ', imoving, ifixed); 
    slicetic = tic;
    movingcloud = pointcell{imoving};
    fixedcloud  = pointcell{ifixed};
    errscurr    = nan(2, 1);
    pcregcurr   = cell(2, 1);

    [tformcurr, pcregcurr{1}, errscurr(1)] = pcregistercpd(movingcloud, fixedcloud, ...
        "Transform","Rigid",'Verbose',false,...
        'Tolerance',1e-6,'OutlierRatio',0.00,'MaxIterations', 100);
    % [tformrev, pcregcurr{2}, errscurr(2)] = pcregistercpd(fixedcloud, movingcloud, ...
    %     "Transform","Rigid",'Verbose',false,...
    %     'Tolerance',1e-6,'OutlierRatio',0.00,'MaxIterations', 100);

    % if errscurr(1)>errscurr(2)
    %     tformcurr = tformrev.invert();
    % end
    tformuse = rigid3dToRigid2d(tformcurr);
    tformarray(islice+1) = affinetform2d(tformuse);

    fprintf('Done! Error match %2.2f. Took %2.2f s\n', errscurr(1), toc(slicetic));
    %--------------------------------------------------------------------------
    % for debugging
    % clf;
    % subplot(1,3,1)
    % plot(movingcloud.Location(:,1), movingcloud.Location(:,2), 'r.',...
    %     fixedcloud.Location(:,1), fixedcloud.Location(:,2), 'k.')
    % axis equal; axis tight; ax = gca; ax.YDir = 'reverse';
    % subplot(1,3,2)
    % plot(pcregcurr{1}.Location(:,1), pcregcurr{1}.Location(:,2), 'r.',...
    %     fixedcloud.Location(:,1), fixedcloud.Location(:,2), 'k.')
    % axis equal; axis tight;ax = gca; ax.YDir = 'reverse';
    % title(sprintf('Error: %3.3f', errscurr(1)))
    % 
    % subplot(1,3,3)
    % plot(movingcloud.Location(:,1), movingcloud.Location(:,2), 'r.',...
    %     pcregcurr{2}.Location(:,1), pcregcurr{2}.Location(:,2), 'k.')
    % axis equal; axis tight;ax = gca; ax.YDir = 'reverse';
    % title(sprintf('Flipped error: %3.3f', errscurr(2)))
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