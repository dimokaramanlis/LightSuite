% Alan Urban Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Gabriel MONTALDO
% Test: Clément BRUNNER, Micheline GRILLET
% May 2020

%% Interpolates and registers a volumetric data with the Allen Mouse Common Coordinate Framework using an affine transformation

% ras=registerData(atlas, x, Transf)
%   atlas, Allen Mouse Common Coordinate Framework provided in the allen_brain_atlas.mat file,
%   x, fus-structure of volume type,
%   Transf, transformation structure obtained with the registering function.
%   ras, the registered anatomy scan

% Example: example03_register_data.m

function ras=registerData(atlas,x,Transf)
Dint=interpolate3D(atlas,x);
T=affine3d(Transf.M);
ref=imref3d(Transf.size);
ras=imwarp(Dint.Data,T,'OutputView',ref);
end




