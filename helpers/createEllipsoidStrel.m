function se = createEllipsoidStrel(semiAxes)
    % createEllipsoidStrel Creates a 3D ellipsoidal structuring element (strel).
    %
    %   se = createEllipsoidStrel(semiAxes) creates a 3D structuring element
    %   object (strel) where the neighborhood is defined by an ellipsoid
    %   centered at the origin with semi-axes lengths given by the input array.
    %
    %   Inputs:
    %       semiAxes - A 3x1 (or 1x3) numeric array [a, b, c] containing the
    %                  semi-axis lengths along the X, Y, and Z dimensions,
    %                  respectively. All values must be >= 0.
    %
    %   Output:
    %       se - A 'strel' object representing the 3D ellipsoid.
    %
    %   Example:
    %       % Create an ellipsoid strel with semi-axes 5, 3, 2
    %       axesVec = [5; 3; 2];
    %       se = createEllipsoidStrel(axesVec);
    %
    %       % Or using a 1x3 array
    %       se_row = createEllipsoidStrel([5, 3, 2]);
    %
    %       % Get the size
    %       disp(size(se.Neighborhood));

    % Validate inputs
    if ~isnumeric(semiAxes) || numel(semiAxes) ~= 3 || any(semiAxes < 0) || ~isvector(semiAxes)
        error('Input must be a 3-element numeric vector (3x1 or 1x3) with non-negative values.');
    end

    % Extract semi-axes
    a = semiAxes(1);
    b = semiAxes(2);
    c = semiAxes(3);

    % Determine the grid size. We need a grid that extends from
    % -ceil(a) to +ceil(a) (and similarly for b and c).
    % The radius for each dimension is ceil(axis)
    rad_x = ceil(a);
    rad_y = ceil(b);
    rad_z = ceil(c);

    % Create coordinate grids centered around 0
    % The length is 2*radius + 1 (to include the 0)
    [X, Y, Z] = meshgrid(-rad_x:rad_x, -rad_y:rad_y, -rad_z:rad_z);

    % Calculate the ellipsoid equation for each point in the grid
    % (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1
    % We add eps (a very small number) to the denominators to avoid
    % division by zero if a, b, or c is 0.
    mask = (X.^2 / (a^2 + eps)) + ...
           (Y.^2 / (b^2 + eps)) + ...
           (Z.^2 / (c^2 + eps)) <= 1;

    % Create the structuring element using the logical mask
    % 'arbitrary' allows us to define a custom neighborhood
    se = strel('arbitrary', mask);
end