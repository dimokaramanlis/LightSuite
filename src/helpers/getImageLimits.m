function imlims = getImageLimits(inputim, alpha)
%GETIMAGELIMITS Summary of this function goes here
%   Detailed explanation goes here

idsuse = find(inputim>0);

minim = quantile(inputim(idsuse),alpha,'all');
maxim = quantile(inputim(idsuse),1-alpha,'all');
maxim  = max(maxim, 5*minim);
imlims = [minim maxim];
end

