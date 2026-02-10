function final_tform = applyAngleToTransform(tformori, Ratlas,cutting_angle_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

saved_vectors = cutting_angle_data.saved_slice_vectors;
valid_vectors = saved_vectors(~cellfun(@(v) any(isnan(v)), saved_vectors));

if isempty(valid_vectors)
    error('cutting_angle_data.mat exists but contains no valid saved vectors.');
end

vector_matrix = cell2mat(valid_vectors);

% --- 2. Calculate rotation in the intuitive GUI coordinate system (AP, ML, DV) ---

% **FIX 1**: Invert the average vector. The GUI camera vector points "back",
% but the slicing normal should point "forward". This fixes the AP inversion.
avg_slice_normal = -mean(vector_matrix, 1);
avg_slice_normal = avg_slice_normal / norm(avg_slice_normal);

original_AP_direction = [1, 0, 0]; % Default AP is along the X-axis in the GUI's space

rotation_axis = cross(original_AP_direction, avg_slice_normal);
rotation_angle_rad = acos(dot(original_AP_direction, avg_slice_normal));

if norm(rotation_axis) > 1e-6 
    k = rotation_axis / norm(rotation_axis);
    K = [  0   -k(3)  k(2); k(3)   0   -k(1); -k(2)  k(1)   0   ];
    R_gui = eye(3) + sin(rotation_angle_rad)*K + (1 - cos(rotation_angle_rad))*(K*K);
else
    R_gui = eye(3) * dot(original_AP_direction, avg_slice_normal);
end

% --- 3. Permute the rotation to match the imwarp world coordinate system ---
P = [0 0 1; 1 0 0; 0 1 0]; % Maps (AP,ML,DV) -> (DV,AP,ML)
R_permuted = P * R_gui * P';

% **FIX 2**: The final matrix needs to be transposed to correct the in-plane
% DV/ML view. This addresses the "transposed view" symptom.
R_final = R_permuted';

% --- 4. Create new transformation preserving the center of rotation ---
original_A = tformori.A;
R_orig = original_A(1:3, 1:3);
T_orig = original_A(1:3, 4);

atlas_center = [mean(Ratlas.XWorldLimits); mean(Ratlas.YWorldLimits); mean(Ratlas.ZWorldLimits)];
T_new = T_orig + (R_orig * atlas_center) - (R_final * atlas_center);

new_A = eye(4);
new_A(1:3, 1:3) = R_final;
new_A(1:3, 4) = T_new;

final_tform = rigidtform3d(new_A);

end