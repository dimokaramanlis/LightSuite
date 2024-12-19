
dp = "D:\staininings\lightsheet\MG391_20230714_redgreen_Mag0.8x_Tile0_Ch488_Sh0_Rot0.tiff";


dz = 5;
dx = 8.2;
dy = 8.2;

cfile = Tiff(dp, 'r');
% get x and y slice numbers;

currvol = tiffreadVolume(dp);
voldown = imresize3(currvol, [1679 1679 974]);
%%
img = double(squeeze(currvol(900,:,:)))';

% vecnoise = sum(img)-movmedian(sum(img),30);
vecnoise = sum(img)-movmedian(sum(img),20);

vecnoise = vecnoise/norm(vecnoise);
imden = img-(img*vecnoise').*vecnoise;

%%
exslice    = squeeze(currvol(:,:,700))';
exslice    = imresize(exslice, 0.333);



%%
% background = single(imgaussfilt(exslice, 100));
% exslice = single(exslice) - background;
[aa,bb,cc] = svd(single(exslice));
Ncomps = 30;
istart = 4;
imden = aa(:,istart:Ncomps)*bb(istart:Ncomps,istart:Ncomps)*cc(:, istart:Ncomps)';
imagesc(imden);
%%
% background = single(medfilt2(exslice, [401 401]));
background = single(imgaussfilt(exslice, 200));
exslice = single(exslice) - background;
imagesc(abs(medfilt2(exslice,[11 11])), [0 1]*30000);
colormap(gray)
axis equal;