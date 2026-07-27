function [slicevol] = generateSliceVolume(sliceinfo, varargin)
%GENERATESLICEVOLUME Assemble a padded, centered slice volume from CZI or TIF files.
%
%   For CZI files: each file may contain multiple scenes (slices).
%   For TIF files: each file in the folder is one slice (no scene selection).
%     Multi-channel TIF files should be multi-page (one page per channel).
%     sliceinfo.pxsize and sliceinfo.channames must be set by the caller.

medfiltwidth = 2*floor((2./sliceinfo.pxsize)/2) + 1;
scalefac     = sliceinfo.pxsize./sliceinfo.px_process;
Nbuff        = ceil(300/sliceinfo.px_process); % 200um for buffer size
size_proc    = ceil(sliceinfo.maxsize.*scalefac);
Nchannels    = numel(sliceinfo.channames);
idx          = 1;

% Detect file format from extension
[~, ~, firstext] = fileparts(sliceinfo.filepaths{1});
is_tif = any(strcmpi(firstext, {'.tif', '.tiff'}));

if nargin > 1
    idreg    = lower(varargin{1});
    regchan  = find(contains(lower(sliceinfo.channames), idreg));
    assert(~isempty(regchan))
else
    regchan      = find(contains(lower(sliceinfo.channames), 'dapi'));
    if isempty(regchan)
        regchan      = find(contains(lower(sliceinfo.channames), 'cy3'));
    end
end

chanids             = 1:Nchannels;
chanids(regchan)    = [];
chanids             = [regchan chanids];
sliceinfo.channames = sliceinfo.channames(chanids);

slicevol     = zeros([size_proc Nchannels sliceinfo.Nslices], 'uint16');
backvalues   = zeros([Nchannels sliceinfo.Nslices], 'uint16');
cropsugg     = zeros(4, sliceinfo.Nslices); % [xmin; xmax; ymin; ymax] per slice
Nfiles       = numel(sliceinfo.filepaths);
maxval = 2^16-1;
%--------------------------------------------------------------------------
dapipx = 2*floor((15./sliceinfo.px_process)/2) + 1;
%--------------------------------------------------------------------------

fprintf('Generating the slice volume by centering slices...\n')
slicetimer   = tic; msg = [];
for ifile = 1:Nfiles

    if is_tif
        % TIF mode: one file = one slice, no scene selection needed
        irel    = [1];
        Nscenes = 1;
    else
        % CZI mode: one file may contain multiple scenes (slices)
        dataim  = BioformatsImage(sliceinfo.filepaths{ifile});
        irel    = sliceinfo.sliceinds{ifile, 2};
        Nscenes = numel(irel);
    end

    for iscene = 1:Nscenes

        if ~is_tif
            dataim.series = irel(iscene);
        end

        %------------------------------------------------------------------
        for icol = 1:Nchannels

            if is_tif
                if Nchannels > 1
                    currim = imread(sliceinfo.filepaths{ifile}, chanids(icol));
                else
                    currim = imread(sliceinfo.filepaths{ifile});
                end
            else
                currim = dataim.getPlane(1, chanids(icol), 1, irel(iscene));
            end

            currim   = medfilt2(currim, medfiltwidth); % to remove salt n' pepper
            currim   = imresize(currim, scalefac(1));
            backval  = quantile(currim(currim>0), 0.01, 'all');
            currim(currim == 0) = backval; % to replace empty tiles

            if sliceinfo.denoisedapi && contains(lower(sliceinfo.channames(icol)), 'dapi')
                currim  = stdfilt(currim, ones(dapipx));
                backval = quantile(currim(currim>0), 0.01, 'all');
                currim(currim == 0) = backval; % to replace empty tiles
            end

            backval  = quantile(currim(currim>0), 0.01, 'all');
            currsize = size(currim);
            padpx    = size_proc - currsize;
            padleft  = floor(padpx/2);

            if icol == 1
                [xrange, yrange] = extractBrainLimits(currim, Nbuff);
                % xrange/yrange are in unpadded-image coords; offset to padded-volume coords
                cropsugg(:, idx) = [xrange(1)   + padleft(2);
                                    xrange(end) + padleft(2);
                                    yrange(1)   + padleft(1);
                                    yrange(end) + padleft(1)];
            end
            currim   = padarray(currim, padleft, backval, 'pre');
            currim   = padarray(currim, padpx - padleft, backval, 'post');
            slicevol(:, :, icol, idx) = currim;
            backvalues(icol, idx)     = backval;

        end
        %------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time/slice %2.2f s. Time elapsed %2.2f s...\n', ...
            idx, sliceinfo.Nslices, toc(slicetimer)/idx, toc(slicetimer));
        fprintf(msg);
        %------------------------------------------------------------------
        idx = idx + 1;
        %------------------------------------------------------------------
    end

end
%--------------------------------------------------------------------------
sliceinfo.size_proc  = size_proc;
sliceinfo.backvalues = backvalues;
%--------------------------------------------------------------------------
% save centering suggestions to decisions file (applied later by alignSliceVolume)
% suggestions are stored in ordering-TIFF pixel coordinates so that
% SliceOrderEditor can display and edit them without knowing the proc resolution
fprintf('Saving centering suggestions to decisions file... ');
ord_scale  = sliceinfo.px_process / sliceinfo.px_register;
ord_size   = ceil(size_proc * ord_scale); % [height, width] of ordering TIFF
cropsugg_ord = max(1, round(cropsugg * ord_scale));
cropsugg_ord([1 2], :) = min(cropsugg_ord([1 2], :), ord_size(2)); % x <= width
cropsugg_ord([3 4], :) = min(cropsugg_ord([3 4], :), ord_size(1)); % y <= height

orderfile = fullfile(sliceinfo.procpath, 'volume_for_ordering_processing_decisions.txt');
if exist(orderfile, 'file')
    T = readtable(orderfile, 'Delimiter', '\t', 'ReadVariableNames', true);
    if height(T) == sliceinfo.Nslices
        T.xmin = cropsugg_ord(1, :)';
        T.xmax = cropsugg_ord(2, :)';
        T.ymin = cropsugg_ord(3, :)';
        T.ymax = cropsugg_ord(4, :)';
    else
        T = makeCropDecisionsTable(sliceinfo.Nslices, cropsugg_ord);
    end
else
    T = makeCropDecisionsTable(sliceinfo.Nslices, cropsugg_ord);
end
writetable(T, orderfile, 'WriteVariableNames', true, 'Delimiter', '\t');
fprintf('Done!\n');
%--------------------------------------------------------------------------
% save volume for processing (centering is NOT yet applied)
fprintf('Saving raw volume... '); tic;
saveLargeSliceVolume(slicevol, sliceinfo.channames, sliceinfo.slicevol);
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------
% save volume for ordering
scalesize   = [ceil(size_proc*sliceinfo.px_process/sliceinfo.px_register) sliceinfo.Nslices];
volproc     = zeros([scalesize(1:2) 3 scalesize(3)], 'uint8');

for ich = 1:min(Nchannels, 3)
    currchan    = squeeze(slicevol(:, :, ich, :));
    currchan    = imresize3(currchan, scalesize);

    backproc    = median(single(backvalues(ich, :)));
    if backproc > 0
        currchan    = (single(currchan) - backproc)./backproc;
    else
        currchan    = single(currchan);
    end
    maxval      = quantile(currchan, 0.99, 'all');
    volproc(:, :, ich, :) = uint8(255*currchan./maxval);
end
options.big      = false;
options.compress = 'lzw';
options.message  = false;
options.color    = true;
if exist(sliceinfo.volorder, 'file')
    delete(sliceinfo.volorder);
end
saveastiff(volproc, sliceinfo.volorder, options);
%--------------------------------------------------------------------------
% let's also save the information file
dpsliceinfo = fullfile(sliceinfo.procpath, 'sliceinfo.mat');
save(dpsliceinfo, 'sliceinfo')
%--------------------------------------------------------------------------
end

function T = makeCropDecisionsTable(Nslices, cropsugg)
T = array2table([(1:Nslices)', zeros(Nslices,1), (1:Nslices)', cropsugg'], ...
    'VariableNames', {'OriginalIndex','FlipState','NewOrderOriginalIndex', ...
    'xmin','xmax','ymin','ymax'});
end
