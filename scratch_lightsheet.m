melastixpath = 'C:\Users\karamanl\Documents\GitHub\matlab_elastix';
yamlmatlabpath = 'C:\Users\karamanl\Documents\GitHub\readyaml';
addpath(genpath(melastixpath), genpath(yamlmatlabpath));



% dpatlas = 'D:\lightsheet\registration\brain_atlas\atlas_template.nii.gz';
dpatlas = 'D:\lightsheet\registration\brain_atlas\gubra_template_olf.nii.gz';
dpano   = 'D:\lightsheet\registration\brain_atlas\gubra_ano_olf.nii.gz';

Vatlas  = niftiread(dpatlas);
Vatlas  = permute(Vatlas, [2 3 1]);
Vatlas  = flip(Vatlas, 3);
Vatlas  = flip(Vatlas, 2);

pxatlas = 20; % in um
%%
dpdata   = 'D:\lightsheet\DK31_Mag0.8x_Tile0_Ch488_Sh1_Rot0.tiff';
dpsave   = 'D:\DATA_folder\Mice\DK031\Anatomy';

dz = 5;
dx = 8.2;
dy = 8.2;
vxsize =  [dx dy dz ];

[voldown, pxcoords] = downsampleLightsheetVolume(dpdata, dpsave, vxsize, 25);

%%
pxcoords = load(fullfile(dpsave, 'downcoords.mat'));
pxcoords = pxcoords.coordsfin;

imfile  = dir(fullfile(dpsave, '*_downsampled.tif'));
voldown = readDownStack(fullfile(imfile.folder, imfile.name));



%%
% Vnew = medfilt3(Vdata, [3 3 5]);
% 
% %%
% 
% islice = 100;
% currslice = Vdata(:,:,islice);
% mask = single(medfilt2(currslice, [101 101]));
% mask = mask > quantile(mask, 0.5,'all');
% subplot(1,2,1)
% imagesc(currslice, [0 5000])
% subplot(1,2,2)
% imagesc(imflatfield(currslice, 100,mask), [0 5000])
%%
% pre-process volume
% imflatfield for gradient correction
% otsu threshold for removing background
% adapthisteq for enhancing contrast
%%
% Vnew     = zeros(size(Vdata),'uint8');
p.Transform='AffineTransform';
p.MaximumNumberOfIterations=400;
p.NumberOfSpatialSamples=300;
p.SP_alpha = 0.2;
tic;
[regfin, stats] = elastix(Vatlas,volnew,[], 'elastix_default.yml','paramstruct',p);
toc;
%%
Rmoving = imref3d(size(Vatlas));
Rfixed  =  imref3d(size(volnew));


[optimizer,metric] = imregconfig("multimodal");
optimizer.MaximumIterations =20;

transtform = imregtform(Vatlas,Rmoving, volnew, Rfixed, ...
    'translation',optimizer, metric,'DisplayOptimization',true); 
[transim, Rtrans]   = imwarp(Vatlas, Rmoving, transtform);
rigidtform = imregtform(transim, Rtrans, volnew, Rfixed, ...
    'rigid',optimizer, metric,'DisplayOptimization',true); 
[rigidim, Rrigid]   = imwarp(transim, Rtrans, rigidtform);


affinetform = imregtform(rigidim, Rrigid, volnew, Rfixed, ...
    'affine',optimizer, metric,'DisplayOptimization',true); 
[affineim, Raffine]   = imwarp(rigidim, Rrigid, affinetform,OutputView=Rfixed);

%%
toshow = imwarp(rigidim, Rrigid, affinetform,OutputView=Rfixed);

viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(tv,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(Vatlas,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);

viewerRegistered2 = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(volnew,Parent=viewerRegistered2,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(toshow,Parent=viewerRegistered2,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);

volshow(volume,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);
