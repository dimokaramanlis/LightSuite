function ptsin = subsampplot(ptsin, npts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
rng(1);
isamp = randperm(size(ptsin,1), min(size(ptsin,1), npts));
ptsin = ptsin(isamp, :);
end