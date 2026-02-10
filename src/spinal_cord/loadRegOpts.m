function regopts = loadRegOpts(dpdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

filecheck = fullfile(dpdata, 'lightsuite/regopts.mat');
if exist(filecheck, 'file')
    regopts = load(filecheck);
else
    error('Could not find regopts in lightsuite folder')
end
end