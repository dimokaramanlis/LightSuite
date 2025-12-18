function tforms = computeStraighteningTransforms(align_out, target_center, target_orientation_deg)
% COMPUTE_STRAIGHTENING_TRANSFORMS
% Creates rigid transforms that straighten the spinal cord.
%
% Logic:
%   1. Translate Slice Center -> Origin
%   2. Rotate (Undo Twist)
%   3. Translate Origin -> Target Center
%
% INPUTS:
%   align_out:     Struct from the GUI (.centers [N x 2], .angles_rad [N x 1])
%   target_center: [x, y] coordinates where the straightened cord should be.
%
% OUTPUT:
%   tforms:        [N x 1] array of rigidtform2d objects.
if nargin < 3
    target_orientation_deg = 90; % Default: Vertical (Top to Bottom)
end


% 1. Retrieve Data
centers = align_out.centers;   
p_ant   = align_out.pred_ant;
p_pos   = align_out.pred_pos;
Nslices = size(centers, 1);

% Handle potential missing data (NaNs)
if any(isnan(centers(:)))
    centers = fillmissing(centers, 'nearest');
    p_ant   = fillmissing(p_ant, 'nearest');
    p_pos   = fillmissing(p_pos, 'nearest');
end

% Pre-allocate
tforms(Nslices, 1) = rigidtform2d; 

% Target Angle in Radians
theta_target = deg2rad(target_orientation_deg);

for i = 1:Nslices
    % --- 1. Current Geometry ---
    cx = centers(i, 1);
    cy = centers(i, 2);
    
    % Calculate Vector: From Anterior -> Posterior
    vec = p_pos(i,:) - p_ant(i,:);
    
    % Calculate Current Angle of this vector
    theta_curr = atan2(vec(2), vec(1));
    
    % Calculate Rotation needed (Target - Current)
    theta_rot = theta_target - theta_curr;
    
    % --- 2. Build Matrices (Pre-Multiplication) ---
    
    % A: Translate Center -> Origin
    T_to_origin = [1,   0,   -cx; 
                   0,   1,   -cy; 
                   0,   0,    1];
     
    % B: Rotate by theta_rot
    % Standard Rotation Matrix for Pre-Multiplication
    % [cos -sin; sin cos]
    c = cos(theta_rot);
    s = sin(theta_rot);
    R = [c,  -s,  0; 
         s,   c,  0; 
         0,   0,  1];
     
    % C: Translate Origin -> Target Center
    tx = target_center(1);
    ty = target_center(2);
    T_to_target = [1,  0,  tx; 
                   0,  1,  ty; 
                   0,  0,  1];
      
    % --- 3. Combine ---
    % v_new = T_target * R * T_origin * v_old
    T_total = T_to_target * R * T_to_origin;
    
    % Create tform object
    tforms(i) = rigidtform2d(T_total);
end