function [registeredMoving, tform] = pcregisterBCPD(moving, fixed, varargin)
%PCREGISTERBCPD Register two point clouds using Bayesian Coherent Point Drift.
%
%   [registeredMoving, transform] = PCREGISTERBCPD(moving, fixed)
%   registers the moving point cloud to the fixed point cloud using the
%   BCPD algorithm with default parameters. It returns the aligned
%   point cloud and a struct containing the transformation parameters.
%
%   [...] = PCREGISTERBCPD(..., Name, Value) specifies additional
%   name-value pair arguments described below.
%
%   Inputs:
%   -------
%   - moving: An M-by-D numeric matrix representing the moving point cloud.
%   - fixed:  An N-by-D numeric matrix representing the fixed point cloud.
%
%   Name-Value Pair Arguments:
%   --------------------------
%   - 'TransformType': Type of transformation.
%     'SimilarityNonrigid' (default) | 'AffineNonrigid' | 'Affine' |
%     'Similarity' | 'Rigid' | 'Nonrigid'
%
%   - 'OutlierRatio': Expected ratio of outliers, w. A value in (0, 1).
%     Default: 0.1
%
%   - 'Lambda': Controls expected deformation length, λ. Smaller is longer.
%     Default: 50.0
%
%   - 'Beta': Controls deformation smoothing range, β.
%     Default: 2.0
%
%   - 'Gamma': Controls rotational invariance, γ. Use ~2-10 for large rotations.
%     Default: 10.0
%
%   - 'ConvergenceTolerance': Convergence tolerance, c.
%     Default: 1e-5
%
%   - 'MaxIterations': Maximum number of iterations, n.
%     Default: 100
%
%   - 'BCPDPath': Full path to the BCPD executable. If empty, it assumes
%     the executable is in the current directory or system path.
%     Default: ''
%
%   - 'Verbose': Display command-line output.
%     true | false (default)
%
%   Outputs:
%   --------
%   - registeredMoving: An M-by-D matrix of the transformed moving point cloud.
%
%   - transform: A struct containing the estimated global transformation.
%     The fields depend on the 'TransformType':
%       - 'Rigid': R (rotation), t (translation)
%       - 'Similarity': s (scale), R (rotation), t (translation)
%       - 'Affine': A (affine matrix), t (translation)
%     For non-rigid types, this is the global part of the transformation.
%
%   Example:
%   --------
%   fixed = rand(100, 3);
%   theta = pi/4;
%   R_true = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
%   moving = fixed * R_true' + [0.5, -0.2, 0.1];
%
%   [registered, transform] = pcregisterBCPD(moving, fixed, ...
%                               'TransformType', 'Rigid', 'Verbose', true);
%
%   disp('Estimated Rotation Matrix:');
%   disp(transform.R);
%   disp('Estimated Translation Vector:');
%   disp(transform.t);

% --- Input Parser ---
p = inputParser;
addRequired(p, 'moving', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'fixed',  @(x) isnumeric(x) && ismatrix(x));

validTransforms = {'similaritynonrigid', 'affinenonrigid', 'affine', 'similarity', 'rigid', 'nonrigid'};
addParameter(p, 'TransformType', 'SimilarityNonrigid', @(x) ismember(lower(x), validTransforms));
addParameter(p, 'OutlierRatio', 0.1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
addParameter(p, 'Lambda', 50.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Beta', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Gamma', 10.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ConvergenceTolerance', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MaxIterations', 500, @(x) isnumeric(x) && isscalar(x) && x > 0 && floor(x)==x);
addParameter(p, 'NormalizeCommon', false, @islogical);
addParameter(p, 'BCPDPath', '', @ischar);
addParameter(p, 'Verbose', false, @islogical);

parse(p, moving, fixed, varargin{:});
params = p.Results;
transform = struct(); % Initialize output transform struct

Ndims = size(moving, 2);
% --- Locate Executable ---
if ~isempty(params.BCPDPath)
    bcpd_exe = params.BCPDPath;
else
    if ispc, bcpd_exe = 'bcpd.exe'; else, bcpd_exe = './bcpd'; end
end
if ~isfile(bcpd_exe) && isempty(dir(which(bcpd_exe)))
    error('BCPD executable not found: %s\nPlease add it to your path or specify ''BCPDPath''.', bcpd_exe);
end

% --- Prepare Temporary Files ---
tmp_prefix = tempname;
if ~exist(tempname, 'dir'), mkdir(tmp_prefix); end
fixedFile    = fullfile(tmp_prefix, 'fixed.txt');
movingFile   = fullfile(tmp_prefix, 'moving.txt');[tmp_prefix, '_moving.txt'];
outputPrefix = fullfile(tmp_prefix, 'output_');

% Define output file paths based on standard BCPD naming
outputFile_y = [outputPrefix, 'y.txt'];
outputFile_s = [outputPrefix, 's.txt'];
outputFile_R = [outputPrefix, 'R.txt'];
outputFile_t = [outputPrefix, 't.txt'];
outputFile_A = [outputPrefix, 'A.txt'];
outputFile_v = [outputPrefix, 'v.txt'];

% Cleanup object to delete all temp files on function exit
cleanupObj = onCleanup(@() rmdir(tmp_prefix, 's'));


writePointsToText(movingFile, moving);
writePointsToText(fixedFile, fixed);

% --- Build BCPD Command String ---
cmd = sprintf('%s -x %s -y %s -o %s', bcpd_exe, fixedFile, movingFile, outputPrefix);

% Map intuitive transform names to BCPD flags
transformMap = containers.Map(...
    {'SimilarityNonrigid', 'AffineNonrigid', 'Affine', 'Similarity', 'Rigid', 'Nonrigid'}, ...
    {'-Tsrn', '-Tan', '-Ta', '-Tsr', '-Tr', '-Tn'});
cmd = [cmd, ' ', transformMap(params.TransformType)];

% Determine which transform files to request
saveFlags = 'cy'; % Always save the resulting point cloud
switch params.TransformType
    case {'Rigid'}
         saveFlags = [saveFlags, 'T'];
    case {'Similarity'}
        saveFlags = [saveFlags, 'T'];
    case {'Affine'}
        saveFlags = [saveFlags, 'T'];
    case {'Nonrigid', 'SimilarityNonrigid', 'AffineNonrigid'}
        saveFlags = [saveFlags, 'Tv'];
end
% cmd = [cmd, '-r1 -p -f0.2 -J300 -Db,5000,0.05 -s', saveFlags];
cmd = [cmd, '-r1 -p -Db,5000,0.05 -s', saveFlags];

% Add other parameters
cmd = [cmd, sprintf(' -w%g -b%g -l%g -g%g', params.OutlierRatio, params.Beta, params.Lambda, params.Gamma)];
cmd = [cmd, sprintf(' -c%g -n%d', params.ConvergenceTolerance, params.MaxIterations)];
if ~params.Verbose, cmd = [cmd, ' -q']; end
if ~params.NormalizeCommon, cmd = [cmd, ' -ux']; end

if params.Verbose
    disp('Executing BCPD command:'); disp(cmd);
end

% --- Execute BCPD ---
[status, cmdout] = system(cmd);
if status ~= 0, error('BCPD execution failed.\nStatus: %d\nOutput:\n%s', status, cmdout); end
if params.Verbose, disp('BCPD Output:'); disp(cmdout); end

% --- Read Results ---
if ~isfile(outputFile_y)
    error('BCPD did not produce the main output file: %s', outputFile_y);
end
registeredMoving = readmatrix(outputFile_y);

% Read transform components if they exist
if isfile(outputFile_s), transform.s = readmatrix(outputFile_s); end
if isfile(outputFile_R), transform.R = readmatrix(outputFile_R); end
if isfile(outputFile_t), transform.t = readmatrix(outputFile_t); end
if isfile(outputFile_A), transform.A = readmatrix(outputFile_A); end
if isfile(outputFile_v), transform.v = readmatrix(outputFile_v); end

%=========================================================================
% fix Rotation matrix
if isfield(transform, 'R')
    % 1. Decompose the CPD estimate using SVD
    [U, ~, V] = svd(transform.R);
    % 2. Reconstruct a strictly orthogonal matrix
    R_clean = U * V';
    % 3. Correct the determinant if it is -1 (reflection)
    % This ensures we have a pure rotation, not a reflection.
    if det(R_clean) < 0
        % Fix the sign of the last column of V (or U)
        V(:,end) = -V(:,end);
        R_clean  = U * V';
    end
    transform.R = R_clean;
end
%==========================================================================

switch params.TransformType
    case {'Rigid'}
        if Ndims == 3
            tform = rigidtform3d(transform.R, transform.t');
        else
            tform = rigidtform2d(transform.R, transform.t');
        end
    case {'Similarity'}
        if Ndims == 3
            tform = simtform3d(transform.s, transform.R, transform.t');
        else
            tform = simtform2d(transform.s, transform.R, transform.t');
        end
    case {'Affine'}
        if Ndims == 3
            tform = affinetform3d([transform.A transform.t]);
        else
            tform = affinetform2d([transform.A transform.t]);
        end
    case {'Nonrigid', 'SimilarityNonrigid'}
        tform.rigid = simtform3d(transform.s, transform.R, transform.t');
        tform.dispv = transform.v;
    case {'AffineNonrigid'}
        tform.affine = affinetform3d([transform.A transform.t]);
        tform.dispv  = transform.v;
end


%==========================================================================
end
