function cf = generateAnnotationMovie(mousestr, foldersave, idim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


savepath      = fullfile('D:\DATA_folder\Mice',mousestr, 'Anatomy');
volume    = readDownStack(fullfile(savepath, 'sample_register_20um.tif'));
transform_params = load(fullfile(savepath, 'transform_params.mat'));


% we permute the volume to match atlas
volume  = permute(volume, transform_params.how_to_perm);
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av     = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
Rmoving  = imref3d(size(av));
Rfixed   = imref3d(size(volume));
avaffine = imwarp(av, Rmoving, transform_params.tform_affine_samp20um_to_atlas_10um_px.invert,...
    'nearest','OutputView',Rfixed);

pathuse = dir(fullfile(savepath, '*atlas_to_samp*.txt'));
pathtrans = fullfile(pathuse.folder, pathuse.name);
anvol = transformAnnotationVolume(pathtrans, avaffine, 20*1e-3);
%==========================================================================
qmax   = single(quantile(volume,0.99,'all'));
qmin   = single(quantile(volume,0.01,'all'))/2;

volume = uint8(255*(single(volume)-qmin)/(qmax-qmin));
sizevol = size(volume);

toplot = true(1, 3);
toplot(idim) = false;
sizeplot = sizevol(toplot);

fcurr = figure('Position',[50 50 sizeplot(2) sizeplot(1)]);


Ny    = size(volume, idim);
ppanel = panel();
ppanel.de.margin = 1;
ppanel.margin = [0 0 0 0];
ppanel.select();ax = gca; ax.Visible = 'off'; axis equal; 
xlim([0.5 sizeplot(2) + 0.5])
ylim([0.5 sizeplot(1) + 0.5])

ax.YDir = 'reverse'; ax.Colormap = gray;

for ii = 1:sizevol(idim)
    ppanel.select(); cla;
    histim  = volumeIdtoImage(volume, [ ii idim]);
    atlasim = volumeIdtoImage(anvol, [ ii idim]);
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim)~=0 & (atlasim > 1);
    [row,col] = ind2sub(size(atlasim), find(av_warp_boundaries));
    image(histim); 
    line(col, row, 'Marker','.','LineStyle','none', 'Color',[1 0.8 0.5],'MarkerSize',4)
    F(ii) = getframe(fcurr) ;
    drawnow;
end

fps = 20;
targetpathcurr = fullfile(foldersave, sprintf('%s_dim%d_registration_video.mp4', mousestr, idim));
writer = VideoWriter(targetpathcurr, 'MPEG-4');
writer.FrameRate = fps;
writer.Quality   = 70;
open(writer)
for i=1:length(F)
    % convert the image to a frame
    writeVideo(writer,  F(i) );
end
%     writeVideo(writer, currvid(:, :, 1, :));
close(writer);
close(fcurr);


end