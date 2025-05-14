function cell_locations = extractCellsFromSliceVolume(opts, ichan)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%------------------------------------------------------------------------
thresuse   = single([0.8 0.5]);
%------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
slicevol = loadLargeSliceVolume(opts.slicevolfin, ichan);
fprintf('Done! Took %2.2f s\n', toc); 
%------------------------------------------------------------------------
[Ny, Nx, Nslices] = size(slicevol);
diaminpx     = opts.celldiam./opts.px_process;
medwithfull  = 2*ceil((2.5*diaminpx)/2) + 1;
matuse       = ones(medwithfull);
Nmed         = floor(sum(matuse, 'all')/2);
%------------------------------------------------------------------------
sigmause   = ceil(diaminpx/4)* [1 1];
cellradius = ceil(diaminpx/2);
%------------------------------------------------------------------------
i0 = 0; % counters
cell_locations = nan(1e6, 5, 'single');
msg = []; tic;
%------------------------------------------------------------------------
for islice = 1:Nslices
    %----------------------------------------------------------------------
    % load slice
    currslice = slicevol(:, :, islice);
    % generate dff image
    backim    = single(ordfilt2(currslice, Nmed, matuse));
    dffim     = gpuArray((single(currslice) - backim)./backim);
    dffim(isnan(dffim)|isinf(dffim)) = 0;
    %----------------------------------------------------------------------
    [ccents, pdiscard] = cellDetector2D(dffim, cellradius, sigmause, thresuse);

    if ~isempty(ccents)
        ikeepx = (ccents(:, 1)> 1) & (ccents(:, 1) < (Nx-1));
        ikeepy = (ccents(:, 2)> 1) & (ccents(:, 2) < (Ny-1));
        ikeep  = ikeepx&ikeepy;
        ccents = ccents(ikeep, :);
        
        if i0+nnz(ikeep)>size(cell_locations,1)
            cell_locations(1e6 + size(cell_locations,1), 1) = 0;
        end
        
        ccents = [ccents(:, 1:2) islice * ones(nnz(ikeep), 1) ccents(:, 3:4)];
        cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
        i0 = i0 + nnz(ikeep);
    end
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d. Points %d. Pdiscard = %2.2f. Time elapsed %2.2f s...\n',...
        islice, Nslices, i0, pdiscard*100, toc);
    fprintf(msg);
    %----------------------------------------------------------------------
end
%----------------------------------------------------------------------
cell_locations = cell_locations(1:i0, :);
% we finally remove weird entries
irem = any(isnan(cell_locations) | isinf(cell_locations), 2);
cell_locations(irem, :) = [];
%--------------------------------------------------------------------------
% after we are done, save cells
if isfield(opts, 'procpath')
    if ~isempty(opts.procpath)
        fsavename = fullfile(opts.procpath, 'cell_locations_sample.mat');
        save(fsavename, 'cell_locations')
    end
end
% delete(opts.fproc)

%--------------------------------------------------------------------------
end