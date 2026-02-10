function [row, col] = getAnnotationBoundaries(annotimg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
annotimg = single(annotimg);
av_warp_boundaries = gradient(annotimg)~=0 & (annotimg > 1);
[row,col] = ind2sub(size(annotimg), find(av_warp_boundaries));
end