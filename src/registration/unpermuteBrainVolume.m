function bvol = unpermuteBrainVolume(tempvol, permvec)
%UNPERMUTEBRAINVOLUME Inverse of permuteBrainVolume.
%   bvol = unpermuteBrainVolume(tempvol, permvec) undoes the flips and
%   permutation applied by permuteBrainVolume(bvol, permvec).

perm_order = abs(permvec);
bvol = tempvol;

% Undo flips first (they were applied after the permute, in permuted space)
for dim = 1:3
    if permvec(dim) < 0
        bvol = flip(bvol, dim);
    end
end

% Undo the permutation
bvol = ipermute(bvol, [perm_order 4]);
end