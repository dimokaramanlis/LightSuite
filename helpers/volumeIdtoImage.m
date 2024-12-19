function  imout = volumeIdtoImage(volume, induse)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S.type = '()';
S.subs = repmat({':'}, 1, 3);  % create a cell array with ':' for all dimensions
S.subs{induse(2)} = induse(1);             % replace the id-th dimension with the index_value
imout = squeeze(subsref(volume, S));



end