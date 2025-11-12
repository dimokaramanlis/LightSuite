function ikeep = selectPtsForAffine(cptshistology)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Npts   = size(cptshistology, 1);
% k      = floor(Npts/2);
% idx    = kmeans(cptshistology, k, 'Replicates', 20);
% ikeep  = accumarray(idx, (1:Npts)', [], @(x) x(1));

% 
% 1. Calculate the centroid of the points
centroid = mean(cptshistology);

% 2. Center the data by subtracting the centroid
centeredPoints = cptshistology - centroid;

% 3. Perform Singular Value Decomposition (SVD)
[~, ~, V] = svd(centeredPoints);

% 4. The direction vector of the line is the first column of U (corresponding to the largest singular value)
lineDirection = V(:, 1);

% Ensure the direction vector points in a consistent direction (e.g., positive z if possible)
if lineDirection(3) < 0
    lineDirection = -lineDirection;
end


% 5. The line origin is the centroid
lineOrigin = centroid'; %' Transpose to get a column vector

% 6. Calculate distances to the line
numPoints = size(cptshistology, 1);
distances = zeros(numPoints, 1);

for i = 1:numPoints
    % Vector from line origin to the point
    vectorToPoint = cptshistology(i, :)' - lineOrigin;  %'

    projection = (lineDirection * (lineDirection' * vectorToPoint));

    % Vector perpendicular to the line (and from the line to the point)
    perpendicularVector = vectorToPoint - projection;

    % Distance is the magnitude of the perpendicular vector
    distances(i) = norm(perpendicularVector);
end

[~, isort] = sort(distances, 'descend');
ikeep      = isort(1:floor(numel(distances)*0.5));
end