function image_out=annotate_volume(image, blob_coords, diameter)
%draws magenta circles on the image

im_double = cast(image, 'double')/255;
spheres = cast(spheres_mask([size(image,1), size(image,2), size(image,3)], blob_coords, diameter), 'double');
color = [1,0,1] ; %magenta
image_out = zeros(size(image,1), size(image,2), size(image,3), 3);
image_out(:,:,:,1)=im_double.*(1-spheres)+color(1)*spheres;
image_out(:,:,:,2)=im_double.*(1-spheres)+color(2)*spheres;
image_out(:,:,:,3)=im_double.*(1-spheres)+color(3)*spheres;
end

function mask = spheres_mask(dims, centers, diameter)
%an array with 1s in the shape of several hollow circles, with centers
% at 'centers' and diamters 'diameters'
mask=zeros(dims);
for i=1:size(centers)
    rad=diameter/2;
    rad_low = (rad)^2;
    rad_high=(rad+1)^2; %outline thickness of 1
    [x,y,z] = ndgrid(-(rad+1):rad+1, -(rad+1):rad+1, -(rad+1):rad+1);
    cx=centers{i,1};
    cy=centers{i,2};
    cz=centers{i,3};
    d=x.^2+y.^2+z.^2;
    try
        mask(cx-(rad+1):cx+rad+1, cy-(rad+1):cy+rad+1, cz-(rad+1):cz+rad+1) = mask( ...
            cx-(rad+1):cx+rad+1, cy-(rad+1):cy+rad+1, cz-(rad+1):cz+rad+1) | (d < rad_high & d>= rad_low);
    catch e
        disp(e.message)
        disp(['probably sphere at ' num2str([cx cy cz]) ' out of bounds'])
    end
end
end
