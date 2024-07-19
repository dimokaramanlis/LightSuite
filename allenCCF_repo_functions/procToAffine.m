function tform = procToAffine(ptransform)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s = ptransform.b;
R = ptransform.T;
T = ptransform.c(1, :) * (1 - s) + ptransform.c(2, :) * s; % Centered translation

% Construct the affine transformation matrix
A = eye(4); % Start with the identity matrix
A(1:3, 1:3) = s * R; % Scale and rotate
A(1:3, 4) = T; % Translate

tform = affine3d(A);
end