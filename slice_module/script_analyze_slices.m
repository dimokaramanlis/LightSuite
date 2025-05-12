% dp            = 'D:\example_charlie\CGF028'; % OG path, good
dp            = 'D:\example_charlie\CGF028';
procpath      = fullfile(dp, 'lightsuite'); makeNewDir(procpath);
dpsliceinfo   = fullfile(procpath, 'sliceinfo.mat'); 
filelistcheck = dir(fullfile(dp, '*.czi'));
filepaths     = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
Nfiles        = numel(filepaths);

px_process     = 2;  % um
px_register    = 20; % um
slicethickness = 150; % um, slice thickness
% here we figure out  sizing
sliceinfo                = getSliceInfo(filepaths);
sliceinfo.procpath       = procpath;
sliceinfo.px_process     = px_process;
sliceinfo.px_register    = px_register;
sliceinfo.slicethickness = slicethickness;
%%
% we first generate the slice volume
[slicevol, sliceinfo] = generateSliceVolume(sliceinfo);
%% reorder slices if needed
InteractiveSliceReorder(fullfile(sliceinfo.procpath, "volume_for_ordering.tiff"))
%% flip slices if needed
SliceFlipper(fullfile(sliceinfo.procpath, "volume_for_ordering.tiff"))

generateReordedVolume(sliceinfo);
%%
InteractiveSliceReorder(fullfile(sliceinfo.procpath, "volume_for_ordering.tiff"))
generateReordedVolume(sliceinfo);
%%
% we align slices
sliceinfo  = load(fullfile(procpath, "sliceinfo.mat"));
sliceinfo  = sliceinfo.sliceinfo;
alignedvol = alignSliceVolume(sliceinfo.slicevol, sliceinfo);
%%

% order slices




%%
clf;
iex = 12;
detailradius = 100; % in um
eximages    = single(gpuArray(slicevol(:,:,regchan,iex+[0 1])));
bval        = median(single(backvalues));
eximages    = (eximages - bval)./bval;
[row, col] = find(spatial_bandpass(eximages(:,:,1), 50, 3, 3, true)>1);
pccurr = pointCloud([gather(col) gather(row) randn(size(col))*32]);
pcdown = pcdownsample(pccurr, 'nonuniformGridSample', 300);

iextest = iex + 1;
[rowtest, coltest] = find(spatial_bandpass(eximages(:,:,2), 50, 3, 3, true)>1);

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

