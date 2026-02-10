function tform = elastix2affine3d(filename)
% ELASTIX2AFFINETFORM3D_PIXEL Converts Elastix transform to MATLAB pixel space
%   tform = elastix2affinetform3d_pixel(filename) reads an Elastix parameter
%   file and returns an affinetform3d object that operates in intrinsic
%   pixel coordinates (1-based).
%
%   It mathematically performs: T_pixel = inv(Grid) * T_physical * Grid

    % 1. Read file content
    txt = fileread(filename);

    % --- HELPER: Parse numeric arrays from text ---
    function val = getParam(name, count)
        pattern = ['\(' name '\s+([^)]+)\)'];
        tokens = regexp(txt, pattern, 'tokens');
        if isempty(tokens)
            error(['Parameter ' name ' not found in file.']);
        end
        str_vals = tokens{1}{1};
        val = sscanf(str_vals, '%f');
        if length(val) ~= count && count ~= -1
            warning(['Expected ' num2str(count) ' values for ' name ', found ' num2str(length(val))]);
        end
    end

    % 2. Extract Physical Transform Parameters (mm)
    % TransformParameters: 12 values (9 rotation + 3 translation)
    params = getParam('TransformParameters', 12);
    
    % CenterOfRotationPoint: 3 values (cx, cy, cz)
    center = getParam('CenterOfRotationPoint', 3);
    
    % Extract Rotation (A) and Translation (t)
    % ITK is row-major, so we reshape and transpose
    A_phys = reshape(params(1:9), [3, 3])'; 
    t_phys = params(10:12);
    
    % Compute the "Effective" Physical Matrix (4x4)
    % T_phys(x) = A(x - c) + t + c  =>  Ax + (t + c - Ac)
    phys_offset = t_phys + center - (A_phys * center);
    
    T_physical = eye(4);
    T_physical(1:3, 1:3) = A_phys;
    T_physical(1:3, 4)   = phys_offset;

    % 3. Extract Grid Parameters (Pixel <-> mm)
    spacing = getParam('Spacing', 3);
    origin  = getParam('Origin', 3);
    
    % Try to find Direction (default to identity if missing)
    try
        direction = getParam('Direction', 9);
        R_dir = reshape(direction, [3, 3])';
    catch
        R_dir = eye(3);
    end

    % 4. Build the Pixel-to-Physical Matrix
    % Formula: x_mm = Origin + Direction * Spacing * (Index_matlab - 1)
    % Note: (Index_matlab - 1) converts 1-based MATLAB to 0-based ITK.
    
    % Scale Matrix (Spacing)
    S_mat = diag(spacing);
    
    % The rotation part of the grid transform
    Grid_Rot = R_dir * S_mat;
    
    % The translation part of the grid transform
    % We must subtract 1 unit of "Grid_Rot" from Origin because MATLAB starts at 1
    Grid_Trans = origin - (Grid_Rot * [1; 1; 1]);
    
    T_pix2phys = eye(4);
    T_pix2phys(1:3, 1:3) = Grid_Rot;
    T_pix2phys(1:3, 4)   = Grid_Trans;

    % 5. Compute Final Composite Transform
    % We want: Pixel_Out = T_total * Pixel_In
    % Step A: Pixel_In -> mm       (T_pix2phys)
    % Step B: Transform in mm      (T_physical)
    % Step C: mm -> Pixel_Out      (inv(T_pix2phys))
    
    T_pixel = T_pix2phys \ (T_physical * T_pix2phys);

    % 6. Create Object
    % affinetform3d expects the premultiply convention, which we maintained.
    tform = affinetform3d(T_pixel);
    
    % --- Debug Output (Optional) ---
    fprintf('Loaded Transform from: %s\n', filename);
    fprintf('Physical Translation (mm): [%.4f, %.4f, %.4f]\n', phys_offset);
    fprintf('Resulting Pixel Translation: [%.4f, %.4f, %.4f]\n', T_pixel(1:3,4));
end