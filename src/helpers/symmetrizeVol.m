function [volsym] = symmetrizeVol(vol, axisuse)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nmid   = size(vol, axisuse)/2;
switch axisuse

    case 1 
        volsym = (vol(1:Nmid, :, :) + flip(vol(:, Nmid+1:end,:, :),1))/2;
    case 2
        volsym = (vol(:, 1:Nmid, :) + flip(vol(:, Nmid+1:end,:),2))/2;
    case 3
        volsym = (vol(:, :, 1:Nmid) + flip(vol(:, :, Nmid+1:end),3))/2;
end
%     
% 
volsym = cat(axisuse, volsym, flip(volsym, axisuse));

end