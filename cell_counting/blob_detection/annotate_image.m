function image_out=annotate_image(image, blob_coords, diameters, varargin)
%draws magenta circles on the image
if size(blob_coords, 2)==3 && all(size(diameters)==1) %2D
    im_double = cast(imadjust(image), 'double')/255;
    diameters=repmat(diameters, [size(blob_coords),1]);
elseif size(blob_coords, 2)==4 %3D
    z=varargin{1};
    im_double = cast(imadjust(image(:,:,z)), 'double')/255;
    d=diameters;
    diameters = zeros([1, size(blob_coords,1)]);
    for i =1:size(blob_coords,1)
    dist = abs(z-blob_coords{i,3}) ;
    if dist>d/2
        diameters(i)=1; %mark for point to indicate out-of-plane sphere
    else
        diameters(i)=round(sqrt((d/2)^2-dist^2))*2;
    end
    end
end
circles = cast(circles_mask([size(image,1), size(image,2)], blob_coords, diameters), 'double');
color = [1,0,1] ; %magenta
image_out = zeros(size(image,1), size(image,2), 3);
image_out(:,:,1)=im_double.*(1-circles)+color(1)*circles;
image_out(:,:,2)=im_double.*(1-circles)+color(2)*circles;
image_out(:,:,3)=im_double.*(1-circles)+color(3)*circles;
end

function mask = circles_mask(dims, centers, diameters)
%an array with 1s in the shape of several hollow circles, with centers
% at 'centers' and diamters 'diameters'
mask=zeros(dims);
for i=1:size(centers)
    if diameters(i)>1
        rad=diameters(i)/2;
        rad_low = (rad)^2;
        rad_high=(rad+1)^2; %outline thickness of 1
        %TODO faster by smaller ndgrid, OR it with relevant part of image only
        [x,y] = ndgrid(-(rad+1):rad+1, -(rad+1):rad+1);
        cx=centers{i,1};
        cy=centers{i,2};
        d=x.^2+y.^2;
        try
        mask(cx-(rad+1):cx+rad+1, cy-(rad+1):cy+rad+1) = mask(cx-(rad+1):cx+rad+1, cy-(rad+1):cy+rad+1) | (d < rad_high & d>= rad_low);
        catch
            disp([num2str(cx-(rad+1)) num2str(cx+rad+1) num2str(cy-(rad+1)) num2str(cy+rad+1)])
        end
    elseif diameters(i)==1 %draw a point
        cx=centers{i,1};
        cy=centers{i,2};
        mask(cx, cy) = 1;
    end
end
end
