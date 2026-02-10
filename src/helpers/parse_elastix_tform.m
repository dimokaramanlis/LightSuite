function tform = parse_elastix_tform(filename)
% READ_ELASTIX_TO_MATLAB Converts Elastix affine file to MATLAB affinetform3d
%
%   tform = read_elastix_to_matlab('affine.txt')
%
%   Returns an affinetform3d object that transforms the Moving Image 
%   to the Fixed Image space (aligning sample to atlas).

    fid = fopen(filename, 'r');
    raw_text = fscanf(fid, '%c');
    fclose(fid);

    % --- 1. Extract Raw Parameters ---
    % Elastix Parameters: [R11 R12 R13 R21 ... Tx Ty Tz]
    params = get_param(raw_text, 'TransformParameters', 12);
    % Center of Rotation: [Cx Cy Cz]
    center = get_param(raw_text, 'CenterOfRotationPoint', 3);
    
    % --- 2. Construct Elastix Matrix (Fixed -> Moving) ---
    % ITK/Elastix uses row-major matrix storage
    R_flat = params(1:9);
    R = reshape(R_flat, [3, 3])'; % Transpose to get standard algebra rotation
    T = params(10:12);
    C = center;
    
    % The Elastix transform is: y = R*(x - C) + T + C
    % We can simplify this to a single 4x4 matrix M where: y = M * x
    % y = R*x - R*C + T + C
    % Effective Translation T_eff = T + C - R*C
    T_eff = T + C - (R * C);
    
    % Construct 4x4 Homogeneous Matrix (Standard Linear Algebra: y = M*x)
    M_elastix = eye(4);
    M_elastix(1:3, 1:3) = R;
    M_elastix(1:3, 4)   = T_eff;
    
    % --- 3. Compute Inverse (Moving -> Fixed) ---
    % Elastix maps Output(Fixed) -> Input(Moving).
    % MATLAB imwarp needs Input(Moving) -> Output(Fixed).
    M_matlab_phys = inv(M_elastix);
    
    % --- 4. Adjust for MATLAB Coordinates (1-based vs 0-based) ---
    % Elastix uses physical coords (origin 0). MATLAB uses intrinsic (origin 1).
    % We must wrap the transform: Shift(-1) -> Transform -> Shift(+1)
    
    shift_to_zero = eye(4);
    shift_to_zero(1:3, 4) = [-1; -1; -1];
    
    shift_to_one = eye(4);
    shift_to_one(1:3, 4) = [1; 1; 1];
    
    % Composite: M_final = Shift(+1) * M_phys * Shift(-1)
    M_final = shift_to_one * M_matlab_phys * shift_to_zero;
    
    % --- 5. Convert to MATLAB Convention (Transposed) ---
    % MATLAB uses post-multiplication: [x y z 1] * M_tform
    % So we transpose our standard algebra matrix.
    tform = affinetform3d(M_final);
    
    % fprintf('Elastix Determinant (Scaling): %.4f\n', det(R));
    % fprintf('MATLAB Inverse Determinant:    %.4f\n', det(M_final(1:3,1:3)));
end

function val = get_param(text, name, count)
    pat = sprintf('\\(%s\\s+([^)]+)\\)', name);
    tokens = regexp(text, pat, 'tokens');
    if isempty(tokens)
        error('Parameter %s not found in file.', name);
    else
        val = sscanf(tokens{1}{1}, '%f');
    end
end