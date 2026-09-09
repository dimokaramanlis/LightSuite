function [XTrain, nrm] = prepareImagesForCNN(cellimages, sigmause, nrm)
% PREPAREIMAGESFORCNN Resizes and stacks cell views for deep learning.
%
% Inputs:
%   cellimages: [Ncells x Nfeatures] flattened raw data
%   sigmause:   [3 x 1] vector of radii
%   nrm:        (optional) [min max] rescaling constants to apply instead of
%               taking them from cellimages. Pass the second output of an
%               earlier call to normalize several batches of cells exactly the
%               same way - the rescaling is global over whatever is handed in,
%               so classifying a large file in chunks needs the constants of
%               the whole file, not of each chunk.
%
% Output:
%   XTrain: [H x W x 3 x Ncells] 4D array for MATLAB CNN training
%   nrm:    [min max] constants actually used


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

if nargin < 3 || isempty(nrm)
    nrm = [min(XTrain, [], 'all'), max(XTrain, [], 'all')];
end
denom = nrm(2) - nrm(1);
if ~isfinite(denom) || denom == 0
    denom = 1;
end
XTrain = (XTrain - nrm(1)) / denom;


end
