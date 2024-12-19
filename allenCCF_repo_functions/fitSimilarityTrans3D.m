
function [tform, mse] = fitSimilarityTrans3D(sourcePoints, targetPoints)
% Ensure the points are in the correct format (Nx3)


% Check input dimensions
if size(sourcePoints, 2) ~= 3 || size(targetPoints, 2) ~= 3
    error('Both sourcePoints and targetPoints must have size Nx3.');
end
if size(sourcePoints, 1) ~= size(targetPoints, 1)
    error('sourcePoints and targetPoints must have the same number of points.');
end

% Compute centroids of each set of points
centroid1 = mean(sourcePoints, 1);
centroid2 = mean(targetPoints, 1);

% Center the points
sourcePoints_centered = sourcePoints - centroid1;
targetPoints_centered = targetPoints - centroid2;

% Compute the cross-dispersion matrix
H = sourcePoints_centered' * targetPoints_centered;

% Perform Singular Value Decomposition (SVD)
[U, ~, V] = svd(H);

% Compute the rotation matrix
R = V * U';
if det(R) < 0
    % Ensure non-reflective rotation
    V(:, end) = -V(:, end);
    R = V * U';
end

% Compute the scale factor
scaleNumerator = sum(sum((targetPoints_centered * R) .* sourcePoints_centered));
scaleDenominator = sum(sum(sourcePoints_centered .^ 2));
s = scaleNumerator / scaleDenominator;

% Compute the translation vector
t = centroid2' - s * R * centroid1';


tform = affinetform3d(simtform3d(s, R, t'));

% assert(size(sourcePoints, 2) == 3, 'Source points must be Nx3');
% assert(size(targetPoints, 2) == 3, 'Target points must be Nx3');
% assert(size(sourcePoints, 1) == size(targetPoints, 1), 'Source and target points must have the same number of rows');
% 
% % Number of points
% n = size(sourcePoints, 1);
% 
% tform = estgeotform3d(sourcePoints,targetPoints,"similarity",'MaxDistance',1000);
% fitgeotrans(sourcePoints, targetPoints)
% % Calculate the mean squared error
transformedSource = transformPointsForward(tform, sourcePoints);
mse = mean(sum((transformedSource - targetPoints).^2, 2));
end