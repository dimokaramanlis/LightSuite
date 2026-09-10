function volout = transformCordImageSlices(cordvol, tforms, raout)
%TRANSFORMCORDIMAGESLICES Warp every slice of a cord volume with its own transform.
%
%   VOLOUT = TRANSFORMCORDIMAGESLICES(CORDVOL, TFORMS, RAOUT) applies TFORMS(k)
%   to slice k of CORDVOL and writes the result into the output frame RAOUT
%   (an imref2d, normally the cross-sectional frame of the atlas). This is how
%   the cord is straightened: the transforms come from
%   COMPUTESTRAIGHTENINGTRANSFORMS, so each section is moved to a common centre
%   and rotated to a common orientation.
%
%   Areas that fall outside the input are filled with the background level of
%   the volume rather than zero, so the straightened volume does not gain hard
%   edges that a later intensity-based registration would try to match.
%
%   See also COMPUTESTRAIGHTENINGTRANSFORMS, CORDSTRAIGHTENPOINTS.

isamps  = randperm(numel(cordvol), min(2e4, numel(cordvol)));
suse    = cordvol(isamps);
backval = mode(suse(suse > 0));

raim   = imref2d(size(cordvol, [1 2]));
volout = zeros([raout.ImageSize numel(tforms)], 'like', cordvol);
for ii = 1:numel(tforms)
    volout(:, :, ii) = imwarp(cordvol(:, :, ii), raim, tforms(ii), ...
        'OutputView', raout, 'FillValues', backval);
end

end
