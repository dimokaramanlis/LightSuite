function [finimage,iy,ix] = blankImage(currimage, blankid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Ny, Nx] =size(currimage);
finimage = 0*currimage;

ix = 1:Nx;
iy = 1:Ny;

if blankid(end)
    ally = linspace(0, Ny/2, 5);
    iy   = round(ally(blankid(1))) + (1:floor(Ny/2));
else
    allx = linspace(0, Nx/2, 5);
    ix = round(allx(blankid(2))) + (1:floor(Nx/2));
end



finimage(iy, ix) = currimage(iy, ix);

end