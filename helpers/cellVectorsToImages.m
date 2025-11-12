function imout = cellVectorsToImages(cellimages, sigmause)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


Ncells    =  size(cellimages, 1);
nk        = [1 3; 2 3; 1 2];
nsigma    = 2*sigmause(nk)+1;
Nperimage = prod(nsigma, 2);
imends   = cumsum(Nperimage);
imstarts = [1; imends(1:2) + 1];

imout = cell(3, 1);

for ii = 1:3
    currim    = cellimages(:, imstarts(ii):imends(ii));
    imout{ii} = reshape(currim, [Ncells nsigma(ii, :)]); 
end
imout = cat(3, imout{:});

end