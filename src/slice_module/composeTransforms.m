function compositeTransform = composeTransforms(transformArray)
% composeTransforms Composes an array of geometric transform objects.
%   compositeTransform = composeTransforms(transformArray)
%   composes the transformations in transformArray in the order they appear.
%   The first transform in the array is applied first to a point/object.
%
%   Input:
%   transformArray: A 1D array or a cell array of transform objects.
%                   Supported objects include rigidtform3d, affinetform3d,
%                   transltform3d, simtform3d, scalartform3d, and their 2D
%                   counterparts (e.g., rigidtform2d). All transforms
%                   must have an '.A' property (transformation matrix) and
%                   preferably a '.Dimensionality' property.
%
%   Output:
%   compositeTransform: A single transform object (e.g., rigidtform3d,
%                       affinetform3d) representing the combined
%                       transformation. The type is chosen based on the
%                       input types and the properties of the resulting matrix.
%
%   Mathematical Convention:
%   If transformArray = {T1, T2, ..., Tn}, and Mi is the matrix for Ti:
%   A point P is transformed as P_final = Mn * ... * M2 * M1 * P.
%   The function computes M_composite = Mn * ... * M2 * M1.
%
%   Example (3D):
%       % Define individual transforms
%       T1_translate = transltform3d([1 2 3]); % Translation
%       S2_scale = scalartform3d(2);           % Uniform scaling by 2
%       R3_rotate = rigidtform3d([0 0 pi/2], [0 0 0]); % Rotation around Z-axis
%
%       % Compose them: T1_translate first, then S2_scale, then R3_rotate
%       transforms3D = {T1_translate, S2_scale, R3_rotate};
%       compTform3D = composeTransforms(transforms3D);
%
%       % Display the resulting matrix
%       disp('Composite 3D Matrix:');
%       disp(compTform3D.A);
%       % Expected type might be affinetform3d due to scaling
%       disp(['Composite 3D Type: ', class(compTform3D)]);
%
%   Example (2D):
%       T1_rigid = rigidtform2d(pi/4, [10 5]); % Rotation and translation
%       A2_affine = affinetform2d([1.5 0.2 0; 0.1 0.5 0; 3 2 1]); % Affine
%
%       transforms2D = {T1_rigid, A2_affine};
%       compTform2D = composeTransforms(transforms2D);
%
%       disp('Composite 2D Matrix:');
%       disp(compTform2D.A);
%       disp(['Composite 2D Type: ', class(compTform2D)]);

    % --- Input Validation and Initialization ---
    if nargin < 1 || isempty(transformArray)
        error('MATLAB:composeTransforms:EmptyInput', 'Input transformArray cannot be empty.');
    end

    % Convert to cell array if it's an object array for uniform access
    if ~iscell(transformArray)
        if isobject(transformArray) && numel(transformArray) > 0
            tempArray = cell(1, numel(transformArray));
            for k_idx = 1:numel(transformArray)
                tempArray{k_idx} = transformArray(k_idx);
            end
            transformArray = tempArray;
        else
            error('MATLAB:composeTransforms:InvalidInputType', ...
                  'Input transformArray must be a cell array of transform objects or an array of transform objects.');
        end
    end

    if ~iscell(transformArray) || isempty(transformArray)
         error('MATLAB:composeTransforms:InvalidInputType', ...
                  'Input transformArray must be a non-empty cell array of transform objects.');
    end

    firstTransform = transformArray{1};
    if ~isobject(firstTransform) || ~isprop(firstTransform, 'A')
         error('MATLAB:composeTransforms:InvalidTransformObject', ...
               'Elements in transformArray must be valid transform objects with an "A" property/method.');
    end

    % Determine dimensionality (2D or 3D)
    dim = 0;
    if isprop(firstTransform, 'Dimensionality')
        dim = firstTransform.Dimensionality;
        if ~(dim == 2 || dim == 3)
            error('MATLAB:composeTransforms:UnsupportedDimensionality', ...
                  'Unsupported dimensionality: %d. Must be 2 or 3.', dim);
        end
    else
        warning('MATLAB:composeTransforms:NoDimProperty', ...
                'First transform object lacks "Dimensionality" property. Inferring from matrix size.');
        try
            matrixA = firstTransform.A; % Access matrix
            matrixSize = size(matrixA);
            if matrixSize(1) == 3 && matrixSize(2) == 3
                dim = 2;
            elseif matrixSize(1) == 4 && matrixSize(2) == 4
                dim = 3;
            else
                error('MATLAB:composeTransforms:UnsupportedMatrixSize', ...
                      'Unsupported matrix size for transformation. Must be 3x3 (2D) or 4x4 (3D).');
            end
        catch ME
            error('MATLAB:composeTransforms:CannotAccessMatrixA', ...
                  'Cannot access .A property of the first transform object: %s', ME.message);
        end
    end
    
    is3D = (dim == 3);
    expectedMatrixSize = dim + 1;

    % Initialize composite matrix with the matrix of the first transform
    try
        compositeMatrix = firstTransform.A;
    catch ME
        error('MATLAB:composeTransforms:CannotAccessMatrixA', ...
              'Failed to get matrix .A from the first transform: %s', ME.message);
    end
    
    if size(compositeMatrix,1) ~= expectedMatrixSize || size(compositeMatrix,2) ~= expectedMatrixSize
        error('MATLAB:composeTransforms:MatrixSizeMismatch', ...
            'Matrix size of first transform (%dx%d) does not match its dimensionality (%d -> %dx%d expected).', ...
            size(compositeMatrix,1), size(compositeMatrix,2), dim, expectedMatrixSize, expectedMatrixSize);
    end


    % --- Determine if all input transforms are strictly rigid types ---
    % (This helps decide if we should attempt to create a rigidtform output)
    allInputTypesWereStrictlyRigid = true; % Assume true initially

    function isStrictlyRigid = checkStrictRigidity(tformObj)
        isStrictlyRigid = false; % Default to false
        if isa(tformObj, 'rigidtform3d') || isa(tformObj, 'rigidtform2d') || ...
           isa(tformObj, 'transltform3d') || isa(tformObj, 'transltform2d')
            isStrictlyRigid = true;
        elseif (isa(tformObj, 'simtform3d') || isa(tformObj, 'simtform2d')) && isprop(tformObj, 'Scale')
            isStrictlyRigid = (abs(tformObj.Scale - 1.0) < 1e-9); % Rigid if scale is 1
        elseif isa(tformObj, 'scalartform3d') && isprop(tformObj, 'Scale')
             isStrictlyRigid = (abs(tformObj.Scale - 1.0) < 1e-9); % Rigid if scale is 1
        % For affinetform, even if it happens to be rigid, we classify it as not "strictly rigid type"
        % to ensure the output is affinetform unless all inputs were explicitly rigid types.
        % If an affinetform that is mathematically rigid is included, the output will be affinetform.
        % This is a conservative choice for type determination.
        end
    end

    if ~checkStrictRigidity(firstTransform)
        allInputTypesWereStrictlyRigid = false;
    end

    % --- Composition Loop ---
    % M_composite = M_n * ... * M_2 * M_1
    % Loop computes: M_new = M_current_transform * M_composite_so_far
    for i = 2:numel(transformArray)
        currentTransform = transformArray{i};

        if ~isobject(currentTransform) || ~isprop(currentTransform, 'A') 
             error('MATLAB:composeTransforms:InvalidTransformObject', ...
                   'Element at index %d is not a valid transform object with an "A" property.', i);
        end

        % Check dimensionality consistency
        if isprop(currentTransform, 'Dimensionality') && currentTransform.Dimensionality ~= dim
            error('MATLAB:composeTransforms:DimensionMismatch', ...
                  'All transforms in the array must have the same Dimensionality. Expected %dD, got %dD at index %d.', ...
                  dim, currentTransform.Dimensionality, i);
        end
        
        try
            nextMatrix = currentTransform.A;
        catch ME
             error('MATLAB:composeTransforms:CannotAccessMatrixA', ...
                   'Failed to get matrix .A from transform at index %d: %s', i, ME.message);
        end

        if size(nextMatrix,1) ~= expectedMatrixSize || size(nextMatrix,2) ~= expectedMatrixSize
             error('MATLAB:composeTransforms:MatrixSizeMismatchLoop', ...
                   ['Matrix size of transform at index %d (%dx%d) does not match expected size %dx%d for %dD transforms.'], ...
                   i, size(nextMatrix,1), size(nextMatrix,2), expectedMatrixSize, expectedMatrixSize, dim);
        end

        compositeMatrix = nextMatrix * compositeMatrix;

        if allInputTypesWereStrictlyRigid && ~checkStrictRigidity(currentTransform)
            allInputTypesWereStrictlyRigid = false;
        end
    end

    % --- Create Output Transform Object ---
    if is3D
        if allInputTypesWereStrictlyRigid
            try
                % Attempt to create a rigidtform3d if all inputs were rigid types
                % The rigidtform3d constructor will validate if the matrix is indeed rigid.
                compositeTransform = rigidtform3d(compositeMatrix);
            catch ME
                warning('MATLAB:composeTransforms:NotRigid3D', ...
                        'Resulting 3D matrix is not strictly rigid or input types included non-rigid. Returning affinetform3d. Original error (if any from rigidtform3d): %s', ME.message);
                compositeTransform = affinetform3d(compositeMatrix);
            end
        else
            % If any input was affine, similarity (scale~=1), scalar (scale~=1), or other, default to affinetform3d
            compositeTransform = affinetform3d(compositeMatrix);
        end
    else % 2D
        if allInputTypesWereStrictlyRigid
            try
                compositeTransform = rigidtform2d(compositeMatrix);
            catch ME
                warning('MATLAB:composeTransforms:NotRigid2D', ...
                        'Resulting 2D matrix is not strictly rigid or input types included non-rigid. Returning affinetform2d. Original error (if any from rigidtform2d): %s', ME.message);
                compositeTransform = affinetform2d(compositeMatrix);
            end
        else
            compositeTransform = affinetform2d(compositeMatrix);
        end
    end

end
