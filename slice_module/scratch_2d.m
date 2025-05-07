dp            = 'D:\example_charlie';
filelistcheck = dir(fullfile(dp, '*.czi'));
filepaths     = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
Nfiles        = numel(filepaths);


px_process  = 2;
px_register = 20;

%%
% here we figure out sizing
sliceinfo    = getSliceInfo(filepaths);
%%
medfiltwidth = 2*ceil((2./sliceinfo.pxsize)/2) + 1;
scalefac     = sliceinfo.pxsize./px_process;
size_proc    = ceil(sliceinfo.maxsize.*scalefac);
Nchannels    = numel(sliceinfo.channames);
slicevol     = zeros([size_proc 4 sliceinfo.Nslices], 'uint16');
idx          = 1;
regchan      = find(strcmp(sliceinfo.channames, 'DAPI'));
Nbuff        = 100;
slicetimer   = tic; msg = [];
for ifile = 1:Nfiles
    dataim = BioformatsImage(filepaths{ifile});
    irel   = sliceinfo.sliceinds{ifile, 2};
    Nscenes = numel(irel);
    for iscene = 1:Nscenes
        dataim.series = irel(iscene);
        currim   = dataim.getPlane(1, regchan, 1, irel(iscene));
        currim   = medfilt2(currim, medfiltwidth);
        currim   = imresize(currim, scalefac(1));
        backval  = quantile(currim, 0.01, 'all');
        %------------------------------------------------------------------
        [xrange, yrange] = extractBrainLimits(currim);
        %------------------------------------------------------------------
        
        yrange = yrange(1)-Nbuff:yrange(2)+Nbuff;
        xrange = xrange(1)-Nbuff:xrange(2)+Nbuff;
        yrange(yrange<1 | yrange > size(currim, 1)) = [];
        xrange(xrange<1 | xrange > size(currim, 2)) = [];

        currim = currim(yrange,xrange);
        currsize = size(currim);
        padpx    = size_proc - currsize;
        padleft  = floor(padpx/2);
        currim   = padarray(currim, padleft, backval, 'pre');
        currim   = padarray(currim, padpx - padleft, backval, 'post');
        
        slicevol(:,:,regchan, idx) = currim;
        % imagesc(currim); drawnow;
        % for icol = 1:Nchannels
        %     currim   = dataim.getPlane(1, icol, 1, irel(iscene));
        % 
        %     % image processing: med filt, resizing, padding
        %     currim   = medfilt2(currim, medfiltwidth);
        %     currim   = imresize(currim, scalefac(1));
        %     currsize = size(currim);
        %     padpx    = size_proc - currsize;
        %     padleft  = floor(padpx/2);
        %     currim   = padarray(currim, padleft, 0, 'pre');
        %     currim   = padarray(currim, padpx - padleft, 0, 'post');
        % 
        %     % save and continue
        %     slicevol(:,:,icol, idx) = currim;
        % end
        %------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d. Time elapsed %2.2f s...\n', ...
            idx, sliceinfo.Nslices, toc(slicetimer));
        fprintf(msg);
        %------------------------------------------------------------------
        idx = idx + 1;
        %------------------------------------------------------------------
    end

end

%%
clf;
iex = 25;
detailradius = 100; % in um

% bb      = corr(single(reshape(slicevol(:,:,4,:), [], sliceinfo.Nslices)));
eximages    = single(gpuArray(slicevol(:,:,4,iex+[0 1])));
[row, col] = find(spatial_bandpass(eximages(:,:,1), 50, 3, 3, true)>1e3);
pccurr = pointCloud([gather(col) gather(row) randn(size(col))*32]);
pcdown = pcdownsample(pccurr, 'nonuniformGridSample', 300);

iextest = iex + 1;
[rowtest, coltest] = find(spatial_bandpass(eximages(:,:,2), 50, 3, 3, true)>1e3);

pccurrtest = pointCloud([gather(coltest) gather(rowtest) randn(size(coltest))*32]);
pcdowntest = pcdownsample(pccurrtest, 'nonuniformGridSample', 300);

coltestflip = size_proc(2) - coltest;
pccurrtestflip = pointCloud([gather(coltestflip) gather(rowtest) randn(size(coltest))*32]);
pcdowntestflip = pcdownsample(pccurrtestflip, 'nonuniformGridSample', 300);

[tform, pcregistered, err] = pcregistercpd(pcdowntest,pcdown, "Transform","Rigid",'Verbose',true,...
    'Tolerance',1e-6,'OutlierRatio',0.00,'MaxIterations', 50);

[tformflip, pcregisteredflip, errflip] = pcregistercpd(pcdowntestflip,pcdown, "Transform","Rigid",'Verbose',true,...
    'Tolerance',1e-6,'OutlierRatio',0.00,'MaxIterations', 50);
fprintf('Normal err = %3.3f, flipped err = %3.3f\n', err, errflip)

subplot(1,3,1)
pcshowpair(pcdowntest, pcdown)
view(2)
subplot(1,3,2)
pcshowpair(pcregistered, pcdown)
view(2)
subplot(1,3,3)
pcshowpair(pcregisteredflip, pcdown)
view(2)
if err<errflip
    tformuse = rigid3dToRigid2d(tform);
    imuse    = eximages(:, :, 2);
else
    tformuse = rigid3dToRigid2d(tformflip);
    imuse    = flip(eximages(:, :, 2), 2);
end

raref        = imref2d(size(imuse));
imregistered = imwarp(imuse, tformuse, 'linear','OutputView',raref);
figure;
imshowpair(eximages(:,:,1), imregistered)

%%
imagesc(eximage,[100 5e4])
hold on;
plot(pcdown.Location(:,1), pcdown.Location(:,2),'r.')

