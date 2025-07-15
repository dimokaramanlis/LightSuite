function [finimage,iy,ix] = blankImage_slice(currimage, blankid, usegauss)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Ny, Nx] = size(currimage);
if usegauss
    sigmause = 2;
    finimage = imgaussfilt(currimage, sigmause);
else
    finimage = 0*currimage;
end

facround = 0.52;
if blankid(1) == 1
    miny = 1;
    maxy = round(Ny*facround);
else
    miny = Ny - round(Ny*facround) + 1;
    maxy = Ny;
end

if blankid(2) == 1
    minx = 1;
    maxx = round(Nx*facround);
else
    minx = Nx - round(Nx*facround) + 1;
    maxx = Nx;
end
iy  = miny:maxy;
ix = minx:maxx;


finimage(iy, ix) = currimage(iy, ix);

end