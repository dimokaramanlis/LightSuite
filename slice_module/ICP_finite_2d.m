function [Points_Moved, M] = ICP_finite_2d(Points_Static, Points_Moving, Options)
% This function ICP_FINITE_2D is a 2D Iterative Closest Point (ICP)
% registration algorithm for point clouds using finite difference methods.
%
% By using finite differences, this function can solve for translation, 
% rotation, resizing, and shearing, unlike standard ICP which is limited
% to rigid transformations.
%
% The function accelerates the registration process by sorting the static 
% points into a grid of overlapping blocks. This allows for a faster search
% for the closest static point to each moving point.
%
% [Points_Moved, M] = ICP_finite_2d(Points_Static, Points_Moving, Options);
%
% Inputs:
%   Points_Static : An N x 2 array with XY points for the registration target.
%   Points_Moving : An M x 2 array with XY points that will be moved and
%                   registered onto the static points.
%   Options       : A struct with registration options:
%     .Registration: 'Rigid'   - Translation and Rotation (default)
%                    'Size'    - Rigid transformations plus resizing
%                    'Affine'  - Translation, Rotation, Resize, and Shear
%     .TolX        : Tolerance for registration position. Default is the
%                    largest side of the bounding box of the points / 1000.
%     .TolP        : Tolerance for the change in distance error. Default is 0.001.
%     .Optimizer   : The optimizer to be used. Options are 'fminlbfgs' (default),
%                    'fminsearch', and 'lsqnonlin'.
%     .Verbose     : If true, displays registration progress. Default is true.
%
% Outputs:
%   Points_Moved  : An M x 2 array with the registered moving points.
%   M             : The 3x3 transformation matrix.
%
% Example:
%   % Create Static Points
%   n_points = 2000;
%   angle_static = linspace(0, 2*pi, n_points);
%   Points_Static = [cos(angle_static)', sin(angle_static)'];
%
%   % Create and Transform Moving Points
%   angle_moving = linspace(0, 2*pi, n_points - 100);
%   Points_Moving = [0.8*cos(angle_moving)', 0.8*sin(angle_moving)'];
%   M_true = [cos(0.2) -sin(0.2) 0.5; sin(0.2) cos(0.2) -0.3; 0 0 1];
%   Points_Moving = [Points_Moving, ones(size(Points_Moving,1),1)] * M_true';
%   Points_Moving = Points_Moving(:,1:2);
%
%   % Register the points
%   Options = struct('Registration','Affine', 'Verbose', true);
%   [Points_Moved, M] = ICP_finite_2d(Points_Static, Points_Moving, Options);
%
%   % Visualize the results
%   figure;
%   hold on;
%   plot(Points_Static(:,1), Points_Static(:,2), 'b.', 'DisplayName', 'Static');
%   plot(Points_Moving(:,1), Points_Moving(:,2), 'r.', 'DisplayName', 'Initial Moving');
%   plot(Points_Moved(:,1), Points_Moved(:,2), 'g.', 'DisplayName', 'Registered');
%   legend;
%   axis equal;
%   title('2D ICP Registration');
%   hold off;

% --- Default Options ---
defaultoptions = struct('Registration','Rigid', 'TolX',0.001, 'TolP',0.001, 'Optimizer','fminlbfgs', 'Verbose', true);
if ~exist('Options','var')
    Options = defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if ~isfield(Options, tags{i})
            Options.(tags{i}) = defaultoptions.(tags{i});
        end
    end
end

% --- Input Validation ---
if size(Points_Static,2) ~= 2
    error('ICP_finite_2d:inputs', 'Points_Static must be an N x 2 matrix.');
end
if size(Points_Moving,2) ~= 2
    error('ICP_finite_2d:inputs', 'Points_Moving must be an M x 2 matrix.');
end

Points_Static = double(Points_Static);
Points_Moving = double(Points_Moving);
Options.Optimizer = lower(Options.Optimizer);

% --- Initial Parameters based on Registration Type ---
switch lower(Options.Registration(1))
    case 'r' % Rigid
        if Options.Verbose, disp('Starting Rigid registration'); end
        scale = [1 1 0.01]; % Translation and Rotation
        par = [0 0 0];
    case 's' % Size
        if Options.Verbose, disp('Starting Size registration'); end
        scale = [1 1 0.01 0.01 0.01]; % Translation, Rotation, and Resize
        par = [0 0 0 100 100];
    case 'a' % Affine
        if Options.Verbose, disp('Starting Affine registration'); end
        scale = [1 1 0.01 0.01 0.01 0.01 0.01]; % Translation, Rotation, Resize, Shear
        par = [0 0 0 100 100 0 0];
    otherwise
        error('ICP_finite_2d:inputs', 'Unknown registration method.');
end

% --- ICP Loop ---
fval_old = inf;
fval_perc = 0;
Points_Moved = Points_Moving;
itt = 0;
maxP = max(Points_Static);
minP = min(Points_Static);
Options.TolX = max(maxP - minP) / 1000;

if Options.Verbose
    disp('    Iteration          Error');
end

% --- Grid-based Point Sorting ---
spacing = size(Points_Static,1)^(1/4);
spacing_dist = max(maxP - minP) / spacing;
xa = minP(1):spacing_dist:maxP(1);
xb = minP(2):spacing_dist:maxP(2);
[x,y] = ndgrid(xa,xb);
Points_Group = [x(:) y(:)];
radius = spacing_dist * sqrt(2);
Cell_Group_Static = cell(1, size(Points_Group,1));

for i=1:size(Points_Group,1)
    mult = 1;
    while isempty(Cell_Group_Static{i})
        check = (Points_Static(:,1) > (Points_Group(i,1) - mult*radius)) & ...
                (Points_Static(:,1) < (Points_Group(i,1) + mult*radius)) & ...
                (Points_Static(:,2) > (Points_Group(i,2) - mult*radius)) & ...
                (Points_Static(:,2) < (Points_Group(i,2) + mult*radius));
        Cell_Group_Static{i} = Points_Static(check,:);
        mult = mult + 1.5;
    end
end

Points_Match = zeros(size(Points_Moved));

while fval_perc < (1 - Options.TolP)
    itt = itt + 1;
    
    for i=1:size(Points_Moved,1)
        Point = Points_Moved(i,:);
        dist = sum((Points_Group - repmat(Point, size(Points_Group,1), 1)).^2, 2);
        [~,j] = min(dist);
        
        Points_Group_Static = Cell_Group_Static{j};
        dist = sum((Points_Group_Static - repmat(Point, size(Points_Group_Static,1), 1)).^2, 2);
        [~,j] = min(dist);
        Points_Match(i,:) = Points_Group_Static(j,:);
    end
    
    % --- Optimization ---
    switch Options.Optimizer
        case 'fminlbfgs'
            optim = struct('Display','off', 'TolX',Options.TolX);
            [par,fval] = fminlbfgs(@(p)affine_registration_error(p, scale, Points_Moving, Points_Match), par, optim);
        case 'fminsearch'
            optim = struct('Display','off', 'TolX',Options.TolX);
            [par,fval] = fminsearch(@(p)affine_registration_error(p, scale, Points_Moving, Points_Match), par, optim);
        case 'lsqnonlin'
            optim = optimset('Display','off', 'TolX',Options.TolX);
            [par,fval] = lsqnonlin(@(p)affine_registration_array(p, scale, Points_Moving, Points_Match), par, [], [], optim);
    end
    
    fval_perc = fval / fval_old;
    if Options.Verbose
        fprintf('     %5.0f       %13.6g\n', itt, fval);
    end
    fval_old = fval;
    
    M = get_transformation_matrix_2d(par, scale);
    Points_Moved = move_points_2d(M, Points_Moving);
end

end

% --- Helper Functions ---

function [e, egrad] = affine_registration_error(par, scale, Points_Moving, Points_Static)
    delta = 1e-8;
    M = get_transformation_matrix_2d(par, scale);
    e = calculate_distance_error_2d(M, Points_Moving, Points_Static);
    if nargout > 1
        egrad = zeros(1, length(par));
        for i = 1:length(par)
            par2 = par;
            par2(i) = par(i) + delta;
            M_delta = get_transformation_matrix_2d(par2, scale);
            egrad(i) = (calculate_distance_error_2d(M_delta, Points_Moving, Points_Static) - e) / delta;
        end
    end
end

function dist_total = calculate_distance_error_2d(M, Points_Moving, Points_Static)
    Points_Moved = move_points_2d(M, Points_Moving);
    dist = sum((Points_Moved - Points_Static).^2, 2);
    dist_total = sum(dist);
end

function e_array = affine_registration_array(par, scale, Points_Moving, Points_Static)
    M = get_transformation_matrix_2d(par, scale);
    Points_Moved = move_points_2d(M, Points_Moving);
    e_array = Points_Moved - Points_Static;
end

function Po = move_points_2d(M, P)
    P_homogeneous = [P, ones(size(P,1), 1)];
    Po_homogeneous = P_homogeneous * M';
    Po = Po_homogeneous(:, 1:2);
end

function M = get_transformation_matrix_2d(par, scale)
    par = par .* scale;
    switch length(par)
        case 3 % Rigid: Translation, Rotation
            M = make_transformation_matrix_2d(par(1:2), par(3));
        case 5 % Size: Translation, Rotation, Resize
            M = make_transformation_matrix_2d(par(1:2), par(3), par(4:5));
        case 7 % Affine: Translation, Rotation, Resize, Shear
            M = make_transformation_matrix_2d(par(1:2), par(3), par(4:5), par(6:7));
    end
end

function M = make_transformation_matrix_2d(t, r, s, h)
    if ~exist('r','var') || isempty(r), r = 0; end
    if ~exist('s','var') || isempty(s), s = [100 100]; end
    if ~exist('h','var') || isempty(h), h = [0 0]; end
    
    s = s / 100; % Scale is given in percentage

    T = [1 0 t(1); 0 1 t(2); 0 0 1];
    R = [cosd(r) -sind(r) 0; sind(r) cosd(r) 0; 0 0 1];
    S = [s(1) 0 0; 0 s(2) 0; 0 0 1];
    H = [1 h(1) 0; h(2) 1 0; 0 0 1];
    
    M = T * R * S * H;
end