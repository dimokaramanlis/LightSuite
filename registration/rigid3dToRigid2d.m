function tform2d_obj = rigid3dToRigid2d(tform3d_obj)
    % rigid3dToRigid2d Converts a rigidtform3d object to a rigidtform2d object.
    %
    % This function extracts the XY translation and the Z-axis rotation (heading)
    % from a 3D rigid transformation to define a 2D rigid transformation.
    %
    % How it works:
    % A rigidtform3d object represents a 3D rigid transformation, defined by a
    % 3x3 rotation matrix (R_3D) and a 3x1 translation vector (t_3D).
    % Its 4x4 homogeneous matrix is T_3D = [R_3D, t_3D; 0 0 0 1].
    %
    % A rigidtform2d object represents a 2D rigid transformation, defined by a
    % 2x2 rotation matrix (R_2D) and a 2x1 translation vector (t_2D).
    % Its 3x3 homogeneous matrix is T_2D = [R_2D, t_2D; 0 0 1].
    %
    % This function determines R_2D and t_2D as follows:
    % 1. t_2D is extracted directly from the X and Y components of t_3D:
    %    t_2D = [t_3D(1); t_3D(2)].
    % 2. R_2D is determined by calculating the effective rotation angle around
    %    the Z-axis (yaw or heading) from R_3D. This angle, theta_z, is
    %    calculated as atan2(R_3D(2,1), R_3D(1,1)).
    %    Then, R_2D = [cos(theta_z) -sin(theta_z);
    %                  sin(theta_z)  cos(theta_z)].
    %
    % This interpretation effectively projects the 3D pose onto the XY plane.
    %
    % Args:
    %   tform3d_obj (rigidtform3d): The input 3D rigid transformation object.
    %
    % Returns:
    %   tform2d_obj (rigidtform2d): The corresponding 2D rigid transformation object.
    %
    % Note: This function requires MATLAB's Navigation Toolbox or a similar toolbox
    %       that defines rigidtform3d and rigidtform2d objects.
    %       The calculation of theta_z using atan2(R_3D(2,1), R_3D(1,1)) is a
    %       standard way to extract yaw assuming a ZYX Euler angle convention
    %       when the pitch angle is not +/- 90 degrees.

    % --- Input Validation ---
    if ~isa(tform3d_obj, 'rigidtform3d')
        error('Input must be a rigidtform3d object.');
    end

    % --- Extract 3D Rotation and Translation ---
    R_3D = tform3d_obj.Rotation;         % 3x3 rotation matrix
    % tform3d_obj.Translation is a 1x3 row vector.
    % We need the X and Y components for the 2D translation.
    translation_3D_row = tform3d_obj.Translation;
    t_2D_col = [translation_3D_row(1); translation_3D_row(2)]; % 2x1 column vector

    % --- Calculate 2D Rotation Angle (theta_z from R_3D) ---
    % This angle represents the rotation in the XY plane (yaw/heading).
    % It's derived from the elements of R_3D that define the orientation
    % of the rotated x-axis in the world XY plane.
    % R_3D = [r11 r12 r13;
    %         r21 r22 r23;
    %         r31 r32 r33]
    % theta_z = atan2(projection of body's x-axis onto world y-axis,
    %                 projection of body's x-axis onto world x-axis)
    % This corresponds to atan2(R_3D(2,1), R_3D(1,1)) for ZYX Euler sequence's Z rotation.
    theta_z = -atan2(R_3D(2,1), R_3D(1,1));

    % --- Construct 2D Rotation Matrix ---
    cos_theta_z = cos(theta_z);
    sin_theta_z = sin(theta_z);
    R_2D = [cos_theta_z, -sin_theta_z;
            sin_theta_z,  cos_theta_z];

    % --- Create rigidtform2d Object ---
    % The rigidtform2d constructor can take R_2D (2x2) and t_2D (2x1 or 1x2)
    % If t_2D_col is 2x1, it should be fine.
    % Alternatively, pass t_2D as a row vector: translation_3D_row(1:2)
    tform2d_obj = rigidtform2d(R_2D, t_2D_col);
    % Or using the row vector for translation:
    % tform2d_obj = rigidtform2d(R_2D, translation_3D_row(1:2));

end