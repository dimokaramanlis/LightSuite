function T = createFlipTransform(volumeSize, flipVector)
    % createFlipTransform Generates a 4x4 transform to flip/rotate points around a volume's center.
    %
    % Inputs:
    %   volumeSize - A 1x3 vector [Ny, Nx, Nz] representing the dimensions
    %                of the 3D volume (MATLAB indexing convention).
    %   flipVector - A 1x3 logical vector [flip_x, flip_y, flip_z] where
    %                true indicates a 180-degree flip/rotation in that axis.
    %
    % Output:
    %   T - A 4x4 homogeneous transformation matrix.
    %
    % Notes:
    %   - The 'flip' is centered on the geometric center of the volume, 
    %     calculated using 1-based indexing (e.g., (N+1)/2).
    %   - A 'Rigid' transform (rotation) occurs when an EVEN number of axes 
    %     are flipped (e.g., [true, true, false] is a 180-deg rotation 
    %     around the Z-axis).
    %   - An 'Improper' transform (reflection) occurs when an ODD number of 
    %     axes are flipped (e.g., [true, false, false]). Your example
    %     is a reflection, not a technically "rigid" transform.
    %
    % Example Usage:
    %   volSize = [100, 120, 80]; % Ny=100, Nx=120, Nz=80
    %
    %   % 1. Identity transform
    %   T_identity = createFlipTransform(volSize, [false, false, false]);
    %
    %   % 2. Flip X-axis only (Reflection, not rigid)
    %   % This matches the user request: [true, false, false] -> x is flipped
    %   T_flip_x = createFlipTransform(volSize, [true, false, false]);
    %
    %   % 3. Flip X and Y axes (180-deg Z-rotation, IS rigid)
    %   T_rot_z = createFlipTransform(volSize, [true, true, false]);
    %

    % 1. Extract dimensions and find center
    % We map volumeSize [Ny, Nx, Nz] to [size_y, size_x, size_z]
    % and then create the center vector in [x, y, z] order.
    size_y = volumeSize(1);
    size_x = volumeSize(2);
    size_z = volumeSize(3);
    
    % Center coordinates (using 1-based indexing)
    center_x = (size_x + 1) / 2;
    center_y = (size_y + 1) / 2;
    center_z = (size_z + 1) / 2;
    
    % Column vector for center [x; y; z]
    center = [center_x; center_y; center_z];

    % 2. Create the scaling/flip matrix (S)
    % Start with identity scale [1, 1, 1]
    scale_factors = [1, 1, 1];
    
    % Set scale to -1 for each axis that needs to be flipped
    % flipVector is [flip_x, flip_y, flip_z]
    if flipVector(1)
        scale_factors(1) = -1; % Flip X
    end
    if flipVector(2)
        scale_factors(2) = -1; % Flip Y
    end
    if flipVector(3)
        scale_factors(3) = -1; % Flip Z
    end
    
    % S is the 3x3 scaling/rotation part of the transform
    S = diag(scale_factors);

    % 3. Calculate the translation part (t)
    % The transform for a point p is: p' = S*p + t
    % We want the flip to be around the center c: p' - c = S * (p - c)
    % p' = S*p - S*c + c
    % p' = S*p + (eye(3) - S) * c
    % So, the translation vector t is (eye(3) - S) * c
    
    t = (eye(3) - S) * center;

    % 4. Assemble the final 4x4 homogeneous matrix
    % T = [ S   t ]
    %     [ 0   1 ]
    T = eye(4);
    T(1:3, 1:3) = S;
    T(1:3, 4) = t;

end