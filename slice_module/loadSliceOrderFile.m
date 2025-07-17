function [outputArg1,outputArg2] = loadSliceOrderFile(sliceinfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

orderfilepath = fullfile(sliceinfo.procpath, 'orderfile.txt');
sliceinfo.ordfilepath = orderfilepath;

if exist(orderfilepath, "file")

end

end