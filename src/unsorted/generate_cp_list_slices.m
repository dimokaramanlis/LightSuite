function cplist = generate_cp_list_slices(volumeuse)
%GENERATE_CP_LIST Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz ] = size(volumeuse);
Nmin          = min([Ny, Nx, Nz]);
assert(Nmin == Ny, 'Slice index is not the first!')

indmat = [1 1; 1 2; 2 1; 2 2];
Ntypes = size(indmat, 1);

% method 3
rng(1);

cplist = [kron((1:Ny)', ones(Ntypes, 1)) kron(ones(Ny,1), indmat)];

cplist = reshape(cplist, Ntypes, Ny, 3);
iperm  = randperm(Ny);
cplist = cplist(:, iperm, :);
cplist = reshape(cplist, [], 3);

% for ii = 1:Ntypes
%     currlist = cplist(ii:Ntypes:end, :);
%     iperm  = randperm(size(currlist, 1));
%     cplist(ii:Ntypes:end, :) = currlist(iperm, :);
% end


end

