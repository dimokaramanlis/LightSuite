function  makeNewDir(dirpath)
%MAKENEWDIR Summary of this function goes here
%   Detailed explanation goes here
if ~exist(dirpath,'dir')
    mkdir(dirpath)
end
end

