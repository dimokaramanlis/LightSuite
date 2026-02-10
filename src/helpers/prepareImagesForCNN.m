function XTrain = prepareImagesForCNN(cellimages, sigmause)
% PREPAREIMAGESFORCNN Resizes and stacks cell views for deep learning.
% 
% Inputs:
%   cellimages: [Ncells x Nfeatures] flattened raw data
%   sigmause:   [3 x 1] vector of radii
%   targetSize: [H W] vector, e.g., [32 32] or [48 48]
%
% Output:
%   XTrain: [H x W x 3 x Ncells] 4D array for MATLAB CNN training


Ncells    =  size(cellimages, 1);
nk        = [1 3; 2 3; 1 2];
nsigma    = 2*sigmause(nk)+1;
Nperimage = prod(nsigma, 2);
imends    = cumsum(Nperimage);
imstarts  = [1; imends(1:2) + 1];
Nmax      = max(nsigma);
XTrain    = zeros([Nmax, 3, Ncells], 'single');

for ii = 1:3
    currim    = cellimages(:, imstarts(ii):imends(ii));
    currshape = nsigma(ii, :);
    currim    = permute(reshape(currim, [Ncells currshape]), [2 3 1]); 
    centval   = floor(Nmax/2) + 1;
    yrange    = centval(1) + (-floor(currshape(1)/2):floor(currshape(1)/2));
    xrange    = centval(2) + (-floor(currshape(2)/2):floor(currshape(2)/2));
    XTrain(yrange, xrange, ii, :) = currim;
end

XTrain = XTrain./sqrt(sum(XTrain.^2, [1 2 3]));
minx   = min(XTrain, [], 'all');
maxx   = max(XTrain, [], 'all');
XTrain = (XTrain - minx) / (maxx - minx);


end