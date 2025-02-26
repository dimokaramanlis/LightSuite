function cplist = generate_cp_list_alt(volumeuse)
%GENERATE_CP_LIST Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nz ] = size(volumeuse);
Nmin          = min([Ny, Nx, Nz]);
minstart      = ceil(Nmin/20);
Nperside      = 20;
Nperm = floor(Nmin*0.9);

indmat = [1 1; 1 2; 2 1; 2 2];
Ntypes = size(indmat, 1);

% method 3
Ncycles = 40;
rng(1);

dataall = cell(3, Ntypes);
for idim = 1:3
    sids = round(linspace(minstart, size(volumeuse, idim)-minstart, Nperside * Ntypes))';
    sids = reshape(sids, Ntypes, Nperside);

    for ii = 1:Ntypes
        dataall{idim, ii} = [sids(ii, :)' idim*ones(Nperside,1) indmat(ii, :).*ones(Nperside,1)];
    end
end

minidx  = min(cellfun(@(x) size(x,1), dataall(:)));
dataall = cellfun(@(x) x(1:minidx, :), dataall, 'UniformOutput',0);

cplist = cat(1, dataall{:});
% iperm  = randperm(size(cplist, 1));
% cplist = cplist(iperm, :);

iperm  = randperm(Nperside);
cplist = reshape(cplist, [Nperside 3 Ntypes Ntypes]);
cplist = permute(cplist(iperm, :, :, :), [4 2 3 1]);
cplist = reshape(cplist, Ntypes, [])';

% method 2
% Ncycles = 10;
% Npercyc = min(floor(Nmin/Ncycles), 15);
% 
% rng(1);
% 
% indskeep = nan(5, Npercyc, 3, Ncycles);
% for icycle = 1:Ncycles
%     Nmove     = (icycle-1)*Npercyc;
%     Nmoveback = Npercyc*Ncycles - Nmove;
%     for idim = 1:3
%         sids = round(linspace(minstart + Nmove, size(volumeuse, idim)-Nmoveback, Npercyc));
%         indskeep(:, :, idim, icycle) = [sids(randperm(Npercyc))' idim * ones(numel(sids),1) randi(5,numel(sids),1)...
%             randi(5,numel(sids),1) rand(numel(sids),1)>0.5]';
%     end
% end
% 
% indskeep = permute(indskeep, [1 3 2 4]);
% cplist = reshape(indskeep, 5, [])';

% method 1
% rng(1);
% for ii = 1:3 %1:3
%     currsize = size(volumeuse, ii);
%     sids = round(linspace(currsize*0.05, currsize*0.95, Nperm));
%     chooselist{ii} = [sids' ii * ones(numel(sids),1) randi(5,numel(sids),1) randi(5,numel(sids),1)...
%         rand(numel(sids),1)>0.5];
% end
% 
% isort = randperm(Nperm);
% cplist = cat(3, chooselist{:});
% cplist = permute(cplist, [2 3 1]);
% cplist = cplist(:, :, isort);
% cplist = reshape(cplist, 5, [])';

end

