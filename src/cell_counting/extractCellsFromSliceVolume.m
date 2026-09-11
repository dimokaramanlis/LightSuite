function cell_locations_all = extractCellsFromSliceVolume(opts, ichan)
%EXTRACTCELLSFROMSLICEVOLUME Detect cells in one or more channels of a slice volume.
%
%   CELL_LOCATIONS = EXTRACTCELLSFROMSLICEVOLUME(OPTS, ICHAN) runs 2-D cell
%   detection on each slice of the aligned volume (opts.slicevolfin) for the
%   channel(s) listed in ICHAN.  ICHAN may be a scalar or a vector.
%
%   Per-channel results are saved to
%       <opts.procpath>/chan<NN>_cell_locations_sample.mat
%
%   Return value:
%     scalar ICHAN  -> Nx5 array  (same behaviour as before)
%     vector  ICHAN -> cell array of Nx5 arrays, one per channel
%------------------------------------------------------------------------
isMultiChan = numel(ichan) > 1;
cell_locations_all = cell(numel(ichan), 1);

if opts.debug
    folderdebug = fullfile(opts.procpath, 'cell_detections');
    makeNewDir(folderdebug)
end
%------------------------------------------------------------------------
for ci = 1:numel(ichan)
    curr_ichan = ichan(ci);
    fprintf('Loading channel %d data... ', curr_ichan); tic;
    slicevol = loadLargeSliceVolume(opts.slicevolfin, curr_ichan);
    if ndims(slicevol) == 4
        slicevol = squeeze(slicevol(:, :, 1, :));
    end
    fprintf('Done! Took %2.2f s\n', toc);
    %----------------------------------------------------------------------
    [Ny, Nx, Nslices] = size(slicevol);
    diaminpx     = opts.celldiam ./ opts.px_process;
    medwithfull  = 2*ceil((2.5*diaminpx)/2) + 1;
    matuse       = ones(medwithfull);
    Nmed         = floor(sum(matuse, 'all') / 2);
    %----------------------------------------------------------------------
    sigmause   = ceil(diaminpx/4) * [1 1];
    cellradius = ceil(diaminpx/2);
    %----------------------------------------------------------------------
    i0 = 0;
    cell_locations = nan(1e6, 5, 'single');
    msg = []; tic;
    %----------------------------------------------------------------------
    for islice = 1:Nslices
        currslice = slicevol(:, :, islice);
        backim    = single(ordfilt2(currslice, Nmed, matuse));
        dffim     = gpuArray((single(currslice) - backim) ./ backim);
        dffim(isnan(dffim) | isinf(dffim)) = 0;
        %------------------------------------------------------------------
        [ccents, pdiscard, imgout] = cellDetector2D(dffim, cellradius, sigmause, opts.thresuse);
        if ~isempty(ccents)
            ikeepx = (ccents(:, 1) > 1) & (ccents(:, 1) < (Nx - 1));
            ikeepy = (ccents(:, 2) > 1) & (ccents(:, 2) < (Ny - 1));
            ikeep  = ikeepx & ikeepy;
            ccents = ccents(ikeep, :);

            if i0 + nnz(ikeep) > size(cell_locations, 1)
                cell_locations(1e6 + size(cell_locations, 1), 1) = 0;
            end

            ccents = [ccents(:, 1:2), islice * ones(nnz(ikeep), 1), ccents(:, 3:4)];
            cell_locations(i0 + (1:nnz(ikeep)), :) = ccents;
            i0 = i0 + nnz(ikeep);
        end
        %------------------------------------------------------------------
        if opts.debug
            pathslice = fullfile(folderdebug, ...
                sprintf('chan%02d_%03d_slice_%d_detections.png', curr_ichan, islice, nnz(ikeep)));
            imtosave = gather(uint8(255 * dffim / opts.thresuse(1)));
            imtosave(imgout) = 255;
            imtosave = cat(3, uint8(imgout*255), imtosave, uint8(imgout*255));
            imtosave = imresize(imtosave, 0.5);
            imwrite(imtosave, pathslice, "png", "BitDepth", 8)
        end
        %------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Chan %d, Slice %d/%d. Points %d. Pdiscard = %2.2f. Time elapsed %2.2f s...\n', ...
            curr_ichan, islice, Nslices, i0, pdiscard*100, toc);
        fprintf(msg);
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    cell_locations = cell_locations(1:i0, :);
    irem = any(isnan(cell_locations) | isinf(cell_locations), 2);
    cell_locations(irem, :) = [];
    %----------------------------------------------------------------------
    if isfield(opts, 'procpath') && ~isempty(opts.procpath)
        fsavename = fullfile(opts.procpath, ...
            sprintf('chan%02d_cell_locations_sample.mat', curr_ichan));
        save(fsavename, 'cell_locations')
    end
    %----------------------------------------------------------------------
    cell_locations_all{ci} = cell_locations;
end
%--------------------------------------------------------------------------
if ~isMultiChan
    cell_locations_all = cell_locations_all{1};
end
%--------------------------------------------------------------------------
end
