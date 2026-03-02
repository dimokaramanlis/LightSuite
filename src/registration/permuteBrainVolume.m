function tempvol = permuteBrainVolume(bvol, permvec)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

perm_order = abs(permvec);
tempvol    = permute(bvol, perm_order);

% Apply flips where value is negative
for dim = 1:3
    if permvec(dim) < 0
        tempvol = flip(tempvol, dim);
    end
end

end