function [slicevol,sliceinfo] = generateSliceVolume(sliceinfo)
%GENERATESLICEVOLUME Summary of this function goes here
%   Detailed explanation goes here

medfiltwidth = 2*floor((2./sliceinfo.pxsize)/2) + 1;
scalefac     = sliceinfo.pxsize./sliceinfo.px_process;
Nbuff        = ceil(200/sliceinfo.px_process); % 200um for buffer size
size_proc    = ceil(sliceinfo.maxsize.*scalefac);
Nchannels    = numel(sliceinfo.channames);
idx          = 1;
regchan      = find(strcmp(sliceinfo.channames, 'DAPI'));
if isempty(regchan)
    regchan      = find(strcmp(sliceinfo.channames, 'Cy3'));
end
chanids             = 1:Nchannels;
chanids(regchan)    = [];
chanids             = [regchan chanids];
sliceinfo.channames = sliceinfo.channames(chanids);

slicevol     = zeros([size_proc Nchannels sliceinfo.Nslices], 'uint16');
backvalues   = zeros([Nchannels sliceinfo.Nslices], 'uint16');
padvalues    = zeros([2 sliceinfo.Nslices]); % for removing later
xrange       = 1:size_proc(2);
yrange       = 1:size_proc(1);
Nfiles       = numel(sliceinfo.filepaths);

fprintf('Generating the slice volume by centering slices...\n')
slicetimer   = tic; msg = [];
for ifile = 1:Nfiles
    dataim = BioformatsImage(sliceinfo.filepaths{ifile});
    irel   = sliceinfo.sliceinds{ifile, 2};
    Nscenes = numel(irel);
    for iscene = 1:Nscenes
        dataim.series = irel(iscene);
        %------------------------------------------------------------------
        for icol = 1:Nchannels
            currim   = dataim.getPlane(1, chanids(icol), 1, irel(iscene));
            currim   = medfilt2(currim, medfiltwidth); % to remove salt n' pepper
            currim   = imresize(currim, scalefac(1));
            backval  = quantile(currim(currim>0), 0.01, 'all');
            if icol == 1
                [xrange, yrange] = extractBrainLimits(currim, Nbuff);
            end

            currim   = currim(yrange,xrange);
            currsize = size(currim);
            padpx    = size_proc - currsize;
            padleft  = floor(padpx/2);
            currim   = padarray(currim, padleft, backval, 'pre');
            currim   = padarray(currim, padpx - padleft, backval, 'post');
            slicevol(:, :, icol, idx) = currim;
            backvalues(icol, idx)     = backval;
        end
        padvalues(:, idx) = padpx;
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
% volume is reduced to save space
sizerem              = floor(min(padvalues, [], 2)/4)*2;
keepx                = (sizerem(2)/2+1):(size_proc(2) - sizerem(2)/2);
keepy                = (sizerem(1)/2+1):(size_proc(1) - sizerem(1)/2);
slicevol             = slicevol(keepy, keepx, :, :);
size_proc            = size_proc - sizerem';
sliceinfo.size_proc  = size_proc;
sliceinfo.backvalues = backvalues;
%--------------------------------------------------------------------------
% save volume for processing
fprintf('Saving volume after centering... '); tic;
saveLargeSliceVolume(slicevol, sliceinfo.channames, sliceinfo.slicevol);
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------
% save volume for ordering
scalesize   = [ceil(size_proc*sliceinfo.px_process/sliceinfo.px_register) sliceinfo.Nslices];
volproc     = zeros([scalesize(1:2) 3 scalesize(3)], 'uint8');
chansmap    = [3 2 1];
for ich = 1:min(Nchannels, 3)
    currchan    = squeeze(slicevol(:, :, ich, :));
    currchan    = imresize3(currchan, scalesize);
    backproc    = reshape(single(backvalues(ich, :)), [1 1 sliceinfo.Nslices]);
    currchan    = (single(currchan) - backproc)./backproc;
    maxval      = quantile(currchan, 0.999, 'all');
    volproc(:, :, chansmap(ich), :) = uint8(255*currchan/maxval);
end
options.big      = false;
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

