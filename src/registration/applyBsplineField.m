function q = applyBsplineField(D, p, direction, varargin)
%APPLYBSPLINEFIELD Move points with a displacement field from READBSPLINEFIELD.
%
%   Q = APPLYBSPLINEFIELD(D, P) returns Q = P + D(P). P is N x 2 or N x 3,
%   1-based voxel coordinates [x y (z)] = [col row (plane)] on the field's
%   grid. D is sampled with linear interpolation; beyond the grid the nearest
%   edge value is used, which only matters for points within half a voxel of
%   the border.
%
%   Q = APPLYBSPLINEFIELD(D, P, 'inverse') solves Q + D(Q) = P, i.e. undoes
%   the forward step exactly, by the fixed-point iteration Q <- P - D(Q). A
%   registration field is smooth next to the voxel size, so this converges in
%   a handful of iterations. Options:
%       'MaxIter' (50)    iteration cap
%       'Tol'     (1e-4)  stop when no point moves more than this (voxels)
%
%   See also READBSPLINEFIELD, BRAINPOINTTRANSFORM.

if nargin < 3 || isempty(direction)
    direction = 'forward';
end
direction = validatestring(direction, {'forward', 'inverse'});

ip = inputParser;
addParameter(ip, 'MaxIter', 50,   @(x) isscalar(x) && x >= 0);
addParameter(ip, 'Tol',     1e-4, @(x) isscalar(x) && x > 0);
parse(ip, varargin{:});

ndim = size(p, 2);
if ndims(D) ~= ndim + 1 || size(D, ndim + 1) ~= ndim
    error('applyBsplineField:dims', ...
        'Points have %d columns but the field has size %s.', ndim, mat2str(size(D)));
end

% The field stays in its own class (single from transformix): griddedInterpolant
% handles single values, and doubling a 20 um brain field would cost ~2 GB for
% no accuracy. Points and results stay double; sampleAll writes into a double
% array, so the iteration below is not limited by single precision.
p = double(p);
F = cell(1, ndim);
for k = 1:ndim
    if ndim == 3
        F{k} = griddedInterpolant(D(:,:,:,k), 'linear', 'nearest');
    else
        F{k} = griddedInterpolant(D(:,:,k), 'linear', 'nearest');
    end
end

switch direction
    case 'forward'
        q = p + sampleAll(F, p);
    case 'inverse'
        q = p;
        for it = 1:ip.Results.MaxIter
            qnew = p - sampleAll(F, q);
            step = max(abs(qnew - q), [], 'all');
            q    = qnew;
            if step < ip.Results.Tol
                break;
            end
        end
        if step >= ip.Results.Tol
            warning('applyBsplineField:noConvergence', ...
                'Inverse did not converge: largest step %.2g voxels after %d iterations.', ...
                step, ip.Results.MaxIter);
        end
end
end

%--------------------------------------------------------------------------
function d = sampleAll(F, p)
d = zeros(size(p));
for k = 1:numel(F)
    if size(p, 2) == 3
        d(:, k) = F{k}(p(:, 1), p(:, 2), p(:, 3));
    else
        d(:, k) = F{k}(p(:, 1), p(:, 2));
    end
end
d(isnan(d)) = 0;
end
