function [tv, av, atlasinfo] = loadCordAtlasVolumes(opts)
%LOADCORDATLASVOLUMES Cord atlas at registration resolution, in sample orientation.
%
%   [TV, AV] = LOADCORDATLASVOLUMES(OPTS) loads the spinal cord template and
%   annotation at OPTS.registres and, if OPTS.tofliprc is true, flips them along
%   the long axis so they run in the same rostrocaudal direction as the sample.
%   Every registration step works in that flipped frame; the flip is undone
%   once, at the very end, by GENERATEREGISTEREDCORDVOLUME and CORDPOINTSTOATLAS.
%
%   The volumes are deliberately not stored inside regopts.mat - they are
%   reloaded here whenever a step needs them, which keeps the options file small
%   and makes the flip state impossible to get out of sync.
%
%   [..., ATLASINFO] = LOADCORDATLASVOLUMES(OPTS) also returns a struct with the
%   region table, the segment table, the native voxel size and the native volume
%   size, i.e. everything needed to relate this grid back to the atlas as it is
%   distributed.
%
%   See also LOADSPINALCORDATLAS, PREPARECORDSAMPLEFORREGISTRATION.

registrationres = opts.registres * [1 1 1];

[tv, av, parcelinfo, segmentinfo, atlasres, nativesize] = ...
    loadSpinalCordAtlas(registrationres);

tofliprc = getOr(opts, 'tofliprc', false);
if tofliprc
    tv = flip(tv, 3);
    av = flip(av, 3);
end

if nargout > 2
    atlasinfo = struct('parcelinfo', parcelinfo, 'segmentinfo', segmentinfo, ...
        'atlasres', atlasres, 'registrationres', registrationres, ...
        'regsize', size(tv), 'nativesize', nativesize, 'tofliprc', tofliprc);
end

end
