function v = symmetrizeVolNan(v, ax)
%SYMMETRIZEVOLNAN Mirror-average a volume about the midline of axis `ax`, NaN-aware.
%   v = symmetrizeVolNan(v, ax) replaces each voxel with the mean of it and its
%   mirror across the midline of dimension `ax` (default 3, the ML axis in atlas
%   space), ignoring NaN. A value present on only one hemisphere is therefore
%   copied to the other (FOV coverage is FILLED, not lost); voxels missing on
%   both sides stay NaN. The result is exactly mirror-symmetric about `ax`, so
%   one hemisphere carries all the (independent) information.
%
%   Unlike SYMMETRIZEVOL, this keeps the volume the SAME size (it does not
%   cat/flip a half back to full) and is NaN-aware, which is what the FOV-gapped
%   fUS volumes need. Valid because the fUSI atlas is mirror-symmetric about the
%   ML midline.
%
%   See also SYMMETRIZEVOL, FUSICROSSMOUSESTRUCTURAL, FUSICOHORTOBJECTPREFERENCE.

if nargin < 2 || isempty(ax), ax = 3; end
nd     = ndims(v);
both   = cat(nd + 1, v, flip(v, ax));
allnan = all(isnan(both), nd + 1);
v      = mean(both, nd + 1, 'omitnan');
v(allnan) = NaN;
end
