function cellimages = getCellImages2D(volumeuse, cinfo, sigmause)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dty = single(-sigmause(1):sigmause(1));
dtx = single(-sigmause(2):sigmause(2));
dtz = single(-sigmause(3):sigmause(3));

[Ny, Nx, Nz] = size(volumeuse);
ccents  =  round(cinfo.WeightedCentroid);
Ncells  = size(ccents, 1);
[xx, yy, zz] = meshgrid(dtx, dty, dtz);

xind   = ccents(:,1)  + xx(:)';
xkeep  = xind>0 & xind <= Nx;
yind   = ccents(:, 2) + yy(:)';
ykeep  = yind>0 & yind <= Ny;
zind   = ccents(:, 3) + zz(:)';
zkeep  = zind>0 & zind <= Nz;

yind(~ykeep) =1;
xind(~xkeep) =1;
zind(~zkeep) =1;

indsout = sub2ind([Ny, Nx, Nz], yind(:), xind(:), zind(:));
Xmat    = reshape(volumeuse(indsout), size(zind));
Xmat(~ykeep | ~xkeep | ~zkeep) = 0;
Xmat = max(reshape(Xmat, Ncells, numel(dty), numel(dtx), numel(dtz)), [], 4);
cellimages = gather(Xmat);


% tic;
% X = zeros(size(cinfo,1), numel(dty), numel(dtx),'single');
% % make the following parallel by indexing simultaneously!!!
% for icell = 1:size(cinfo,1)
%     ccent  =  round(cinfo.WeightedCentroid(icell,:));
%     xind = ccent(1) + dtx;
%     yind = ccent(2) + dty;
%     zind = ccent(3) + dtz;
%     zind(zind<1 | zind > size(volumeuse, 3)) = [];
%     ix = xind>0 & xind <= size(volumeuse, 2);
%     iy = yind>0 & yind <= size(volumeuse, 1);
%     X(icell,iy, ix) = max(volumeuse(yind(iy), xind(ix), zind), [], 3);
% end
% cellimages = reshape(X, size(cinfo,1), numel(dty)*numel(dtx));
% toc;


end