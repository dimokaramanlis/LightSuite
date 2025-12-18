function volout = tranformCordImagesSlices(cordvol, tforms, raout)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

isamps  = randperm(numel(cordvol), min(2e4, numel(cordvol)));
suse    = cordvol(isamps);
backval = mode(suse(suse > 0));

raim   = imref2d(size(cordvol, [1 2]));
volout = zeros([raout.ImageSize numel(tforms)], "uint16");
for ii = 1:numel(tforms)
    volout(:, :, ii)  = imwarp(cordvol(:, :, ii), raim, tforms(ii), ...
        'OutputView',raout, 'FillValues', backval);
end



end