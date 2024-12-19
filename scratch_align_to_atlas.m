
% dpimage = 'S:\ElboustaniLab\#SHARE\Data\AM097 (Npx4)\Anatomy\Axiocam Fluo\20231118\Slide3_Slice2.czi';
dpimage = 'S:\ElboustaniLab\#SHARE\Data\AM090 (Npx3)\Anatomy\AxioSlideScanner\2023_11_24__0001_AM090_1.czi';

dataim  = BioformatsImage(dpimage);
regim   = dataim.getPlane(1, 2, 1, 1);
dpatlas = 'D:\AllenAtlas\template_volume_10um.npy';
atlas   = readNPY(dpatlas);



%%
tv = readNPY(dpatlas); % grey-scale "background signal intensity"
av = readNPY('D:\AllenAtlas\annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
st = loadStructureTree('D:\AllenAtlas\structure_tree_safe_2017.csv'); % a table of what all the labels mean

file_save_location = 'D:\AllenAtlas'; % where will the probe locations be saved
probe_name = 'test'; % name probe to avoid overwriting
f = figure;
f = allenAtlasBrowser(f, tv, av, st, file_save_location, probe_name,'transverse');
%%
%Get rough mask
f1 = figure('Position',[50 50 1200 800]);
imagesc(regim, getImageLimits(regim, 0.01)); colormap(gray)
axis equal; axis tight;
ax = gca; ax.Visible = 'off'; 
ax.Title.Visible = 'on';
title('Use freehand drawing to encircle the slice')
ptsdraw = drawfreehand(ax);
imask = createMask(ptsdraw);
close(f1);
%%
maskedregim = single(regim).*imask;

[irow, icol] = find(imask);
yimlims = [min(irow) max(irow)];
ximlims = [min(icol) max(icol)];

maskedregim = maskedregim(yimlims(1):yimlims(2), ximlims(1):ximlims(2));
% maskedregim = flip(maskedregim, 2);
imagesc(maskedregim, getImageLimits(maskedregim, 0.01));colormap(gray)
axis equal; axis tight;
ax = gca; ax.Visible = 'off'; 
ax.Title.Visible = 'on';
%%

% Nfilt    = floor(30/dataim.pxSize(1)/2)*2+1;
imfiltered  = imsharpen(maskedregim,'Radius',40,'Amount',0.8,'Threshold',0);
% imfiltered  = imgaussfilt(maskedregim, [Nfilt]);
%%
subplot(1,2,1)
imagesc(imfiltered, getImageLimits(imfiltered, 0.01));colormap(gray)

axis equal; axis tight;
ax = gca; ax.Visible = 'off'; 
ax.Title.Visible = 'on';
%%
subplot(1,2,2)
imagesc(atlas(:,:,580)');colormap(gray)

axis equal; axis tight;
ax = gca; ax.Visible = 'off'; 
ax.Title.Visible = 'on';
% irem = imfiltered> quantile(imfiltered,0.999,'all');
% imfiltered(irem) = median(imfiltered, 'all');
%%
limsalign = getImageLimits(imfiltered, 0.01);

imtoalign = (imfiltered-limsalign(1))/range(limsalign);
atlastemplate = atlas(:,:,580)';
limsatlas =  single(getImageLimits(atlastemplate, 0.01));

imatlas  = (single(atlastemplate) - limsatlas(1))/range(limsatlas);
%%
[selectedMovingPoints, selectedFixedPoints] = cpselect(imtoalign, imatlas ,'Wait',true);
%%
RA = imref2d(size(imatlas));
RB = imref2d(size(imtoalign));

tform = fitgeotrans(selectedMovingPoints, selectedFixedPoints, 'affine');
[imregistered,RB]= imwarp(imtoalign, RB, tform);
%%
% boost contrast

minim = quantile(imfiltered,0.01,'all');
maxim = quantile(imfiltered,0.99,'all');
regim = regim - quantile(regim,0.01,'all')