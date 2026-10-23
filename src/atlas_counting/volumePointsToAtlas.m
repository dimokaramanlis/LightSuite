function finalpts = volumePointsToAtlas(inputpts, trstruct)
%VOLUMEPOINTSTOATLAS Map raw sample points to 10 um atlas voxels.
%
%   FINALPTS = VOLUMEPOINTSTOATLAS(INPUTPTS, TRSTRUCT) is kept for older
%   scripts and is the same as
%       brainPointTransform(INPUTPTS, TRSTRUCT, 'raw2atlas')
%   Columns 1:3 of INPUTPTS are raw [x y z] pixels; the rest are carried over.
%
%   See also BRAINPOINTTRANSFORM, TRANSFORMPOINTSTOATLAS.

finalpts = brainPointTransform(inputpts, trstruct, 'raw2atlas');
end
