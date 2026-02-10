function [tform, mse] = fitRigidTrans3D(sourcePoints, targetPoints)
%fitRigidTrans3D Finds the rigid transformation to fit source points to target points.
%   [tform, mse] = fitRigidTrans3D(sourcePoints, targetPoints) calculates the
%   rigid transformation (rotation and translation) that best aligns the
%   sourcePoints to the targetPoints in a least-squares sense.
%
%   Inputs:
%   sourcePoints - An Nx3 matrix of 3D coordinates.
%   targetPoints - An Nx3 matrix of 3D coordinates, corresponding to sourcePoints.
%
%   Outputs:
%   tform - A rigidtform3d object representing the optimal transformation.
%   mse   - The mean squared error between the transformed source points
%           and the target points.
%
%   This implementation is based on the Kabsch algorithm.

% Step 1: Validate Inputs
% Ensure the points are in the correct format (Nx3) and have the same number of points.
assert(size(sourcePoints, 2) == 3, 'Source points must be an Nx3 matrix.');
assert(size(targetPoints, 2) == 3, 'Target points must be an Nx3 matrix.');
assert(size(sourcePoints, 1) == size(targetPoints, 1), 'Source and target point sets must have the same number of points.');

% Step 2: Calculate Centroids
% The centroid is the mean position of all the points in the set.
centroidSource = mean(sourcePoints, 1);
centroidTarget = mean(targetPoints, 1);

% Step 3: Center the Point Sets
% Subtract the centroid from each point set to move them to the origin.
% This separates the translation component from the rotation component.
centeredSource = sourcePoints - centroidSource;
centeredTarget = targetPoints - centroidTarget;

% Step 4: Calculate the Covariance Matrix
% The covariance matrix H captures the correlation between the centered point sets.
% It is calculated as H = P_source_centered' * P_target_centered.
H = centeredSource' * centeredTarget;

% Step 5: Compute the Optimal Rotation using Singular Value Decomposition (SVD)
% Decompose the covariance matrix H into U, S, and V.
% H = U * S * V'
% The optimal rotation matrix R is given by R = V * U'.
[U, ~, V] = svd(H);
R = V * U';

% Step 6: Handle Special Case (Reflection)
% If the determinant of the rotation matrix is -1, it's a reflection (an improper rotation).
% This can happen in the presence of noise or symmetric data. To ensure a
% proper rotation, we flip the sign of the last column of V and recalculate R.
if det(R) < 0
    %fprintf('Reflection detected, correcting...\n');
    V(:, 3) = V(:, 3) * -1;
    R = V * U';
end

% Step 7: Calculate the Optimal Translation
% The translation vector 't' is the difference between the target centroid
% and the rotated source centroid.
% t = centroid_target' - R * centroid_source'
t = centroidTarget' - R * centroidSource';

% Step 8: Create the rigidtform3d Transformation Object
% Use MATLAB's built-in object to represent the rigid transformation.
% The rigidtform3d object stores the rotation (R) and translation (t).
% The constructor expects the translation as a 1x3 row vector.
tform = rigidtform3d(R, t');

% Step 9: Calculate the Mean Squared Error (MSE)
% Apply the found transformation to the original source points.
transformedSourcePoints = (R * sourcePoints' + t)';
% The MSE is the average of the squared Euclidean distances between
% the transformed source points and the original target points.
mse = mean(sum((transformedSourcePoints - targetPoints).^2, 2));

end
