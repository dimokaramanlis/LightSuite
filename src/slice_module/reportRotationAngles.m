function [rx,ry,rz, tstr] = reportRotationAngles(rotationMatrix, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 2. Extract the 3x3 rotation matrix from the rigidtform3d object.

% 3. Manually convert the rotation matrix to Euler angles ('ZYX' convention).
% This section replaces the rotm2eul function.
% Based on the rotation matrix R:
%   R(3,1) = -sin(rotY)
%   R(3,2) = cos(rotY)*sin(rotX)
%   R(3,3) = cos(rotY)*cos(rotX)
%   R(2,1) = cos(rotY)*sin(rotZ)
%   R(1,1) = cos(rotY)*cos(rotZ)

if nargin<2
    verbose = true;
else
    verbose = varargin{1};
end
% Check for gimbal lock: This occurs when rotY is +/- 90 degrees,
% which makes cos(rotY) = 0, causing division by zero in standard formulas.
% This corresponds to R(3,1) being -1 or 1.
if abs(rotationMatrix(3,1)) > 0.99999
    % Gimbal Lock Case
    rotY_rad = -asin(rotationMatrix(3,1));
    rotZ_rad = 0; % Conventionally, set Z rotation (yaw) to 0
    rotX_rad = atan2(-rotationMatrix(1,2), -rotationMatrix(1,3)); % Calculate X rotation
else
    % Standard Case
    rotY_rad = -asin(rotationMatrix(3,1));
    % Use atan2 for numerical stability
    rotX_rad = atan2(rotationMatrix(3,2), rotationMatrix(3,3));
    rotZ_rad = atan2(rotationMatrix(2,1), rotationMatrix(1,1));
end

% The resulting angles correspond to sequential rotations around Z, Y, and X axes.
euler_angles_rad = [rotZ_rad, rotY_rad, rotX_rad];


% 4. Convert the angles from radians to degrees.
euler_angles_deg = rad2deg(euler_angles_rad);


% 5. Assign the results to variables corresponding to your defined axes.
% Based on the 'ZYX' sequence:
% - The first angle is rotation about the Z-axis.
% - The second angle is rotation about the Y-axis.
% - The third angle is rotation about the X-axis.
rz    = euler_angles_deg(1);
ry    = euler_angles_deg(2);
rx    = euler_angles_deg(3);


% 6. Display the final results in degrees.
tstr = sprintf('Global rotation around ML: %.1f deg, DV: %.1f deg, AP: %.1f deg', rz,rx,ry);
if verbose
    fprintf('%s\n', tstr);
end

end