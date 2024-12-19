function [tform, mse] = fitAffineTrans3D(sourcePoints, targetPoints)
% Ensure the points are in the correct format (Nx3)
assert(size(sourcePoints, 2) == 3, 'Source points must be Nx3');
assert(size(targetPoints, 2) == 3, 'Target points must be Nx3');
assert(size(sourcePoints, 1) == size(targetPoints, 1), 'Source and target points must have the same number of rows');

% Number of points
n = size(sourcePoints, 1);

% Augment the source points with a column of ones to account for translation
augmentedSource = [sourcePoints, ones(n, 1)];

% Set up the least-squares problem: AX = B
% A is the augmented source points, X is the affine transformation matrix, B is the target points
A = augmentedSource;
B = targetPoints;

% Solve for X using the backslash operator (A\B)
X = A \ B;

% Construct the affine transformation matrix
affineMatrix = eye(4);
affineMatrix(1:3, :) = X';

% Create the affine3d object
tform = affinetform3d(affineMatrix);
mse   = mean(sum((A*X - B).^2,2));
end