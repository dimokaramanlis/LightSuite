function tform = parse_elastix_tform_all(filename)
% PARSE_ELASTIX_TFORM  Convert any *linear* Elastix transform to affinetform3d.
%
%   tform = parse_elastix_tform('TransformParameters.0.txt')
%
%   Supports every linear Elastix transform:
%       TranslationTransform  (3 params)
%       EulerTransform        (6 params)  - rigid
%       SimilarityTransform   (7 params)  - rigid + isotropic scale
%       AffineTransform       (12 params)
%   Non-linear transforms (BSpline, spline-kernel, ...) cannot be represented
%   as an affinetform3d and raise an error.
%
%   Returns an affinetform3d that maps the MOVING image into the FIXED image
%   space (aligning sample to atlas), in MATLAB *intrinsic* pixel coordinates,
%   ready for  imwarp(movingVol, tform, 'OutputView', imref3d(size(fixedVol))).
%
%   All Elastix linear transforms share the same physical-space model
%       y = A*(x - C) + T + C
%   and only differ in how the 3x3 matrix A is built from the parameters.
%
%   COORDINATE HANDLING
%   The transform is defined in Elastix *physical/world* coordinates. This
%   function reads Spacing / Origin / Direction from the file (they describe
%   the FIXED grid) and converts intrinsic <-> physical properly, so it is
%   correct for non-unit spacing (your file uses Spacing = 0.05).
%
%   ASSUMPTIONS (edit here if they don't hold for your data)
%     * The MOVING image is sampled on the same grid geometry (Spacing,
%       Origin, Direction) as the FIXED image. That geometry is not stored in
%       the transform file; when both images share a resolution (the usual
%       atlas<->sample case) this is exact. If they differ, pass the moving
%       geometry in where S is built below.
%     * Elastix axis order (i,j,k) maps directly to MATLAB (x,y,z). If you had
%       to permute your volumes for the affine version, keep doing the same.

    txt = fileread(filename);

    % --- 1. Read fields (with sensible defaults for optional ones) ----------
    ttype      = get_str(txt, 'Transform');                 % e.g. 'EulerTransform'
    p          = get_vec(txt, 'TransformParameters');       % all optimisable params
    C          = get_vec_or(txt, 'CenterOfRotationPoint', [0;0;0]);
    spacing    = get_vec_or(txt, 'Spacing',   [1;1;1]);
    origin     = get_vec_or(txt, 'Origin',    [0;0;0]);
    dirFlat    = get_vec_or(txt, 'Direction', [1;0;0;0;1;0;0;0;1]);
    Direction  = reshape(dirFlat, [3,3])';                  % row-major -> matrix
    computeZYX = strcmpi(get_str_or(txt, 'ComputeZYX', 'false'), 'true');

    % --- 2. Build the 3x3 linear part A and translation T (physical space) --
    [A, T] = build_linear(ttype, p, computeZYX);

    % --- 3. Elastix forward map: FIXED(phys) -> MOVING(phys) ---------------
    %     y = A*(x - C) + T + C  =>  T_eff = T + C - A*C
    T_eff             = T + C - A*C;
    M_fixed_to_moving = [A, T_eff; 0 0 0 1];

    % --- 4. Invert: imwarp needs MOVING(phys) -> FIXED(phys) ---------------
    M_moving_to_fixed = M_fixed_to_moving \ eye(4);   % == inv(M_fixed_to_moving)

    % --- 5. Physical <-> pixel conversion ---------------------------------
    %   ITK/Elastix:  phys = Origin + Direction*diag(Spacing) * index0based
    %   MATLAB intrinsic is 1-based, so subtract 1 to reach the 0-based index.
    S = eye(4);
    S(1:3,1:3) = Direction * diag(spacing);   % (moving grid, assumed == fixed grid)
    S(1:3,4)   = origin;

    shift_to_zero = eye(4);  shift_to_zero(1:3,4) = [-1; -1; -1];  % 1-based -> 0-based
    shift_to_one  = eye(4);  shift_to_one(1:3,4)  = [ 1;  1;  1];  % 0-based -> 1-based

    % intrinsic(moving) -> 0idx -> phys -> transform -> phys -> 0idx -> intrinsic(fixed)
    M_final = shift_to_one / S * M_moving_to_fixed * S * shift_to_zero;

    % --- 6. affinetform3d uses the premultiply convention  y = M*[x;y;z;1] --
    %     so M_final is passed directly (no transpose).
    tform = affinetform3d(M_final);
end

% =======================================================================
function [A, T] = build_linear(ttype, p, computeZYX)
% Build the 3x3 linear part A and 3x1 translation T for a linear transform.
    t = lower(ttype);
    if contains(t,'spline') || contains(t,'kernel') || contains(t,'weightedcombination')
        error('parse_elastix_tform:nonlinear', ...
            'Transform "%s" is non-linear and cannot be an affinetform3d.', ttype);

    elseif contains(t,'translation')
        need(p, 3, ttype);
        A = eye(3);
        T = p(1:3);

    elseif contains(t,'euler')                     % [ax ay az tx ty tz]
        need(p, 6, ttype);
        A = rot_euler(p(1), p(2), p(3), computeZYX);
        T = p(4:6);

    elseif contains(t,'similarity')                % [q1 q2 q3 tx ty tz s]
        need(p, 7, ttype);
        A = p(7) * rot_versor(p(1), p(2), p(3));   % isotropic scale * rotation
        T = p(4:6);

    elseif contains(t,'affine')                    % [a11..a33(row-major) tx ty tz]
        need(p, 12, ttype);
        A = reshape(p(1:9), [3,3])';               % row-major -> matrix
        T = p(10:12);

    else
        error('parse_elastix_tform:unknown', ...
            'Unrecognised transform type "%s".', ttype);
    end
    T = T(:);
end

function R = rot_euler(ax, ay, az, computeZYX)
% ITK/Elastix Euler3D. Default order is ZXY (ComputeZYX "false"): R = Rz*Rx*Ry.
    Rx = [1 0 0;  0 cos(ax) -sin(ax);  0 sin(ax) cos(ax)];
    Ry = [cos(ay) 0 sin(ay);  0 1 0;  -sin(ay) 0 cos(ay)];
    Rz = [cos(az) -sin(az) 0;  sin(az) cos(az) 0;  0 0 1];
    if computeZYX
        R = Rz*Ry*Rx;
    else
        R = Rz*Rx*Ry;                              % elastix/ITK default
    end
end

function R = rot_versor(x, y, z)
% ITK versor (unit-quaternion vector part) -> rotation matrix. w = scalar part.
    w = sqrt(max(0, 1 - x^2 - y^2 - z^2));
    R = [1-2*(y^2+z^2),  2*(x*y-w*z),    2*(x*z+w*y);
         2*(x*y+w*z),    1-2*(x^2+z^2),  2*(y*z-w*x);
         2*(x*z-w*y),    2*(y*z+w*x),    1-2*(x^2+y^2)];
end

function need(p, n, ttype)
    if numel(p) < n
        error('parse_elastix_tform:params', ...
            '%s expects %d parameters but found %d.', ttype, n, numel(p));
    end
end

% ---- field readers ----------------------------------------------------
function v = get_vec(txt, name)
    tok = regexp(txt, ['\(' name '\s+([^)]*)\)'], 'tokens', 'once');
    if isempty(tok), error('Parameter %s not found in file.', name); end
    v = sscanf(tok{1}, '%f');
end

function v = get_vec_or(txt, name, default)
    tok = regexp(txt, ['\(' name '\s+([^)]*)\)'], 'tokens', 'once');
    if isempty(tok), v = default(:); else, v = sscanf(tok{1}, '%f'); end
end

function s = get_str(txt, name)
    tok = regexp(txt, ['\(' name '\s+"([^"]*)"\)'], 'tokens', 'once');
    if isempty(tok), error('Parameter %s not found in file.', name); end
    s = tok{1};
end

function s = get_str_or(txt, name, default)
    tok = regexp(txt, ['\(' name '\s+"([^"]*)"\)'], 'tokens', 'once');
    if isempty(tok), s = default; else, s = tok{1}; end
end
