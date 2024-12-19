function cplist = generate_cp_list_adv(volumeuse)
%GENERATE_CP_LIST Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz ] = size(volumeuse);
Nmin = min([Ny, Nx, Nz]);
Nperm =  floor(Nmin*0.9);


chooselist = cell(3,1);
rng(1);
for ii = 1:3 %1:3
    currsize = size(volumeuse, ii);
    sids = round(linspace(currsize*0.1, currsize*0.9, Nperm));
    chooselist{ii} = [sids' ii * ones(numel(sids),1)];
end

isort = randperm(Nperm);
cplist = cat(3, chooselist{:});
cplist = permute(cplist, [2 3 1]);
cplist = cplist(:, :, isort);
cplist = reshape(cplist, 2, [])';


sids  = round(linspace(Ny*0.05, Ny*0.95, 10));
cptop = [sids' ones(numel(sids),1)];
cplist = cat(1, cptop, cplist);

end

