function [axis, angle] = rotm2axang_custom(R)
    % Custom implementation to convert a 3D rotation matrix to axis-angle
    % representation without requiring any toolboxes.
    % R: 3x3 rotation matrix
    % axis: 1x3 unit vector representing the axis of rotation
    % angle: scalar angle of rotation in radians

    % Calculate the angle of rotation
    angle = acos((trace(R) - 1) / 2);

    % If the angle is 0, the axis is undefined.
    % We can return any unit vector, e.g., [1 0 0] and an angle of 0.
    if abs(angle) < 1e-6
        axis = [1 0 0];
        angle = 0;
        return;
    end

    % If the angle is pi (180 degrees), sin(angle) is 0, so we have a
    % different way to calculate the axis.
    if abs(angle - pi) < 1e-6
        % Find the axis from R + I = 2*u*u'
        [v, ~] = eigs(R + eye(3), 1); % The axis is the eigenvector for eigenvalue 1
        axis = v';
        angle = pi;
        return;
    end
    
    % The general case
    rx = R(3,2) - R(2,3);
    ry = R(1,3) - R(3,1);
    rz = R(2,1) - R(1,2);

    axis = (1 / (2 * sin(angle))) * [rx, ry, rz];
end