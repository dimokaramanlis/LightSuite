%{
    MATLAB Prototype: Point Cloud Registration via Mutual Information Maximization
    
    This script implements the registration method we discussed.
    It finds a rigid transformation (6-DOF) to align a "moving" cloud (X)
    to a "fixed" cloud (Y) by maximizing the mutual information (MI)
    between their local density features.
    
    Dependencies:
    - Deep Learning Toolbox (for dlarray, dlgradient, dlfeval, adamupdate)
    - Statistics and Machine Learning Toolbox (for pdist2, pcshow)
    - A CUDA-enabled GPU is highly recommended.
%}

clear; close all; clc;

disp("Starting Mutual Information Registration Prototype...");

%% 1. Setup: Create Synthetic Point Clouds
% We'll create a box 'Y' and a transformed copy 'X'
N_y = 2000; % Number of points in fixed cloud
N_x = 1800; % Number of points in moving cloud

% Create a "box" point cloud for Y
Y = rand(N_y, 3, 'single') - 0.5;
Y_edges = Y(abs(Y(:,1)) > 0.4 | abs(Y(:,2)) > 0.4 | abs(Y(:,3)) > 0.4, :);
Y = Y_edges;
N_y = size(Y, 1);
fprintf("Generated fixed cloud 'Y' with %d points.\n", N_y);

% Define a ground-truth transformation (rodriguez_vec, translation_vec)
true_rvec = [0.1, -0.2, 0.15];
true_tvec = [0.05, -0.05, 0.1];
true_R = rotationVectorToMatrix(true_rvec);
true_t = single(true_tvec);

% Create the moving cloud 'X' by transforming Y and adding noise
% --- FIX: Ensure N_x is not greater than the new N_y after filtering ---
N_x = min(N_x, N_y); 
if N_x == 0
    error("Filtering 'Y' resulted in 0 points. Cannot create 'X'. Adjust filtering criteria or increase N_y.");
end
% --- End Fix ---
X = Y(1:N_x, :) * true_R' + true_t;
X = X + (rand(N_x, 3, 'single') - 0.5) * 0.02; % Add small noise
fprintf("Generated moving cloud 'X' with %d points.\n", N_x);

% --- Move data to GPU ---
try
    g = gpuDevice;
    fprintf("GPU detected: %s\n", g.Name);
    X_gpu = gpuArray(X);
    Y_gpu = gpuArray(Y);
    useGPU = true;
catch
    warning("GPU not detected or Parallel Computing Toolbox not available. Running on CPU.");
    X_gpu = X;
    Y_gpu = Y;
    useGPU = false;
end

% Visualize the initial state (optional, requires Statistics Toolbox)
if license('test', 'Statistics_Toolbox')
    figure;
    pcshow(Y, 'green', 'MarkerSize', 50);
    hold on;
    pcshow(X, 'red', 'MarkerSize', 50);
    title('Initial Alignment');
    legend('Fixed Cloud (Y)', 'Moving Cloud (X)', 'Location', 'NorthWest');
    axis equal;
end

%% 2. Setup: Hyperparameters and Query Points

% --- Feature Hyperparameters ---
% This is the 'sigma' for the spatial KDE (your "histogram distance")
% It defines the scale of the local neighborhood.
param.feature_sigma = 0.1;

% --- MI Estimator Hyperparameters ---
% This is the bandwidth for the 2D MI KDE.
% It smooths the (g_Y, g_X) feature distribution.
param.mi_bandwidth = 0.05;

% --- Query Points ---
% We define the features at these locations.
% Using a random subset of Y is a good, simple strategy.
K = 8000; % Total number of query points
query_indices = randperm(N_y, min(N_y, K));
P_gpu = Y_gpu(query_indices, :);
numQueryPoints = size(P_gpu, 1); % <-- Get the actual number of query points
fprintf("Using %d query points (subset of Y).\n", numQueryPoints);

% --- Optimization Hyperparameters ---
batchSize = 1024; % Number of query points per gradient step
numIterations = 200;
learnRate = 0.005; % Adam learn rate

% --- FIX: Ensure batchSize is not larger than the number of query points ---
if batchSize > numQueryPoints
    warning('Hyperparameters: batchSize (%d) is larger than numQueryPoints (%d). Truncating batchSize to %d.', batchSize, numQueryPoints, numQueryPoints);
    batchSize = numQueryPoints;
end
% --- End Fix ---

%% 3. Precomputation: Fixed Cloud Features

% We can pre-compute the density features for the fixed cloud 'Y'
% at all query points 'P'. This doesn't change during optimization.
fprintf("Pre-computing density features for fixed cloud 'Y'...\n");
% --- FIX: Ensure inputs to computeDensityFeatures are dlarray ---
% The helper function expects dlarray so it can use extractdata
P_dl_precompute = dlarray(P_gpu);
Y_dl_precompute = dlarray(Y_gpu);
g_Y_dl = computeDensityFeatures(P_dl_precompute, Y_dl_precompute, param.feature_sigma);
% --- End Fix ---
disp("Precomputation complete.");

%% 4. Optimization: Find Transformation

% --- Initialize Transformation Parameters ---
% 6-DOF: [rx, ry, rz, tx, ty, tz]
% We use 'dlarray' to track gradients.
T_params = zeros(6, 1, 'single');
T_params_dl = dlarray(T_params);
if useGPU
    T_params_dl = gpuArray(T_params_dl);
end

% --- Initialize Adam Optimizer States ---
avgG = [];
% ... (Adam optimizer setup) ...

% --- Prepare 'dlarray' versions of data for the loop ---
% We pass all of 'X' each time.
X_dl = dlarray(X_gpu); % --- FIX: Remove dimension labels ('CB') ---

% Pre-allocate for plotting loss
lossHistory = zeros(numIterations, 1);

disp("Starting optimization loop...");
tic;
for i = 1:numIterations
    % 1. Sample a minibatch of query points
    [P_batch_dl, g_Y_batch_dl] = sampleMinibatch(P_gpu, g_Y_dl, batchSize, useGPU); % <-- FIX: Pass g_Y_dl
    
    % 2. Evaluate the objective function and compute gradients
    %    We use dlfeval to run the function 'objectiveFunction'
    %    and automatically get the gradient of its output (loss)
    %    w.r.t. its input (T_params_dl).
    [loss, grads] = dlfeval(@objectiveFunction, T_params_dl, X_dl, P_batch_dl, g_Y_batch_dl, param);
    
    % 3. Update parameters using Adam
    [T_params_dl, avgG, avgGSq] = adamupdate(T_params_dl, grads, avgG, avgGSq, i, learnRate);
    
    % 4. Log and Display Progress
    lossHistory(i) = gather(extractdata(loss));
    if mod(i, 20) == 0 || i == 1
        fprintf('Iter: %03d | Loss (Neg. MI): %f | rvec: [%.3f, %.3f, %.3f] | tvec: [%.3f, %.3f, %.3f]\n', ...
            i, lossHistory(i), ...
            gather(extractdata(T_params_dl(1:3))), ...
            gather(extractdata(T_params_dl(4:6))));
    end
end
toc;
disp("Optimization finished.");

%% 5. Results

% --- Get final transformation ---
final_T_params = gather(extractdata(T_params_dl));
final_rvec = final_T_params(1:3)';
final_tvec = final_T_params(4:6)';
final_R = rotationVectorToMatrix(final_rvec);

fprintf("\n--- Results ---\n");
fprintf("True R-Vec:     [%.4f, %.4f, %.4f]\n", true_rvec);
fprintf("Estimated R-Vec:  [%.4f, %.4f, %.4f]\n\n", final_rvec);
fprintf("True T-Vec:     [%.4f, %.4f, %.4f]\n", true_tvec);
fprintf("Estimated T-Vec:  [%.4f, %.4f, %.4f]\n", final_tvec);

% --- Apply final transformation ---
X_transformed = X * final_R' + final_tvec;

% --- Visualize final result ---
if license('test', 'Statistics_Toolbox')
    figure;
    pcshow(Y, 'green', 'MarkerSize', 50);
    hold on;
    pcshow(X_transformed, 'magenta', 'MarkerSize', 50);
    pcshow(X, 'red', 'MarkerSize', 20, 'Marker', '.');
    title('Final Alignment');
    legend('Fixed Cloud (Y)', 'Transformed Cloud (X_{final})', 'Original Cloud (X_{init})', 'Location', 'NorthWest');
    axis equal;
    
    figure;
    plot(lossHistory);
    title('Loss Curve (-MI)');
    xlabel('Iteration');
    ylabel('Loss');
    grid on;
end

%% 6. Local Helper Functions
% These functions are defined locally so the script is self-contained.
% They are written to be compatible with dlarray.

function loss = objectiveFunction(T_params_dl, X_dl, P_batch_dl, g_Y_batch_dl, param)
% --- FIX: This function should only return the loss. ---
% dlfeval automatically computes and returns the gradients externally.
%
% This is the function that dlgradient will differentiate.
% It computes the (negative) MI from the transformation parameters.
    
    % 1. Apply current transformation T to the moving cloud X
    %    T(X) = X*R' + t
    T_X_dl = applyTransform(X_dl, T_params_dl);
    
    % 2. Compute density features for the *transformed* moving cloud, g_X
    %    at the query point batch locations P_batch.
    g_X_batch_dl = computeDensityFeatures(P_batch_dl, T_X_dl, param.feature_sigma);
    
    % 3. Estimate Mutual Information from the (g_Y, g_X) pairs
    %    We want to MAXIMIZE MI, so our "loss" to MINIMIZE is -MI.
    loss = -estimateMI(g_Y_batch_dl, g_X_batch_dl, param.mi_bandwidth);
    
    % Gradients are computed automatically by dlgradient.
% ... (objectiveFunction definition) ...
end

function [P_batch_dl, g_Y_batch_dl] = sampleMinibatch(P_gpu, g_Y_dl, batchSize, useGPU) % <-- FIX: Argument name
% Selects a random minibatch of query points and their precomputed features.
    
    numQueryPoints = size(P_gpu, 1);
    % We no longer need the min() check here because we fixed batchSize above
    idx = randperm(numQueryPoints, batchSize); 
    
    P_batch = P_gpu(idx, :);
    
    % --- FIX: g_Y_dl is already a dlarray, so indexing gives a dlarray ---
    g_Y_batch_dl = g_Y_dl(idx, :); 
    
    % Convert P_batch to dlarray for this batch (P_gpu was numeric)
    P_batch_dl = dlarray(P_batch);
    % --- End Fix ---
    
    % No need for 'if useGPU', dlarray will be on GPU if P_batch was.
end


function T_X_dl = applyTransform(X_dl, T_params_dl)
% Applies a 6-DOF transform (as dlarray) to a set of points (as dlarray).
    
    % Extract parameters
    rvec = T_params_dl(1:3); % Rotation vector
    tvec = T_params_dl(4:6); % Translation vector
    
    % --- FIX: Use a custom, differentiable Rodrigues' formula ---
    % Convert rotation vector to 3x3 matrix (differentiable)
    % The built-in rotationVectorToMatrix does not support dlarray.
    R = differentiableRotationVectorToMatrix(rvec);
    
    % Apply transform: T(X) = X * R + t
    % We use R directly (not R') because our ground truth generation:
    % true_R = rotationVectorToMatrix(rvec) (which is R_vision')
    % X = Y * true_R' + t (which is Y * (R_vision')' + t = Y * R_vision + t)
    % Our new function returns R_vision, so we multiply by R.
    % --- End Fix ---
    
    % We need to reshape tvec from [3,1] to [1,3] for broadcasting
    tvec_row = permute(tvec, [2, 1]); % [1, 3]
    
    % X_dl is (N_x, 3), R is (3, 3)
    T_X_dl = X_dl * R + tvec_row;
end

function g = computeDensityFeatures(P_dl, X_dl, sigma)
% Computes the "soft-count" density features:
% g(p) = sum_j( exp(-||p_i - x_j||^2 / (2*sigma^2)) )
%
% P_dl: Kx3 query points
% X_dl: Nx3 cloud points
% g:    Kx1 density features
    
    % --- FIX: pdist2 does not support dlarray directly ---
    % Extract underlying data (e.g., gpuArray) for pdist2
    P_data = extractdata(P_dl);
    X_data = extractdata(X_dl);
    
    % pdist2 is compatible with gpuArray, but not dlarray
    % D2(i, j) = ||p_i - x_j||^2
    D2_data = pdist2(P_data, X_data, 'squaredeuclidean');
    
    % Promote the numeric result back to a dlarray for gradient tracking
    D2 = dlarray(D2_data);
    % --- End Fix ---
    
    % Gaussian kernel
    sigma2 = 2 * sigma^2;
    G = exp(-D2 / sigma2);
    
    % Sum over all cloud points 'x_j' to get the density at each 'p_i'
    g = sum(G, 2); % Result is Kx1
end

function mi = estimateMI(U_dl, V_dl, h)
% Estimates Mutual Information from two Kx1 vectors of features.
% MI(U, V) = H(U) + H(V) - H(U, V)
% using a simple Parzen-window (KDE) entropy estimator.
%
% U_dl: Kx1 (e.g., g_Y)
% V_dl: Kx1 (e.g., g_X)
% h:    KDE bandwidth
    
    K = size(U_dl, 1);
    h2 = 2 * h^2;
    
    % --- FIX: pdist2 does not support dlarray directly ---
    % Extract underlying data for pdist2
    U_data = extractdata(U_dl);
    V_data = extractdata(V_dl);

    % --- H(U) ---
    D2_U_data = pdist2(U_data, U_data, 'squaredeuclidean');
    D2_U = dlarray(D2_U_data); % Promote back to dlarray
    p_u = sum(exp(-D2_U / h2), 2); % Kx1
    % Add epsilon for numerical stability before log
    H_U = -mean(log(p_u + 1e-9)); 
    
    % --- H(V) ---
    D2_V_data = pdist2(V_data, V_data, 'squaredeuclidean');
    D2_V = dlarray(D2_V_data); % Promote back to dlarray
    p_v = sum(exp(-D2_V / h2), 2); % Kx1
    H_V = -mean(log(p_v + 1e-9));
    
    % --- H(U, V) ---
    UV_dl = [U_dl, V_dl]; % Kx2
    UV_data = extractdata(UV_dl); % Extract the combined Kx2 data
    D2_UV_data = pdist2(UV_data, UV_data, 'squaredeuclidean');
    D2_UV = dlarray(D2_UV_data); % Promote back to dlarray
    % --- End Fix ---
    
    p_uv = sum(exp(-D2_UV / h2), 2); % Kx1. Note: h is used for both dims.
    H_UV = -mean(log(p_uv + 1e-9));
    
    % We ignore the constant log(K*h...) terms as they
    % cancel out or don't affect the gradient w.r.t. T_params.
    mi = H_U + H_V - H_UV;
end

% --- NEW FUNCTION ---
function R = differentiableRotationVectorToMatrix(rvec)
% Implements Rodrigues' formula for a 3x1 dlarray rotation vector
% to produce a 3x3 differentiable rotation matrix.
% This is the standard formula R = I + sin(theta)*K + (1-cos(theta))*K^2
% which corresponds to p_new = p_old * R
    
    % Add a small epsilon for numerical stability in norm and division
    epsilon = 1e-12; 
    
    % Calculate angle (theta)
    theta_sq = sum(rvec.^2);
    theta = sqrt(theta_sq);
    
    % Calculate axis (k)
    % We use (theta + epsilon) to avoid NaN gradients at theta=0
    k = rvec / (theta + epsilon); 
    
    kx = k(1);
    ky = k(2);
    kz = k(3);
    
    % --- FIX: Assign numeric data, not dlarray, to gpuArray ---
    rvec_data = extractdata(rvec); 
    
    % Extract the underlying numeric data from the dlarray components
    kx_data = extractdata(kx);
    ky_data = extractdata(ky);
    kz_data = extractdata(kz);
    
    % Skew-symmetric matrix K (as a numeric gpuArray first)
    K = zeros(3, 3, 'like', rvec_data); 
    K(1,2) = -kz_data;
    K(1,3) = ky_data;
    K(2,1) = kz_data;
    K(2,3) = -kx_data;
    K(3,1) = -ky_data;
    K(3,2) = kx_data;
    
    % Promote K to a dlarray to track gradients from here
    K_dl = dlarray(K);
    
    % K^2 (now a dlarray operation)
    K2 = K_dl * K_dl; 
    
    % Rodrigues' formula: R = I + sin(theta)*K_dl + (1-cos(theta))*K^2
    I = eye(3, 'like', rvec_data); 
    % --- End Fix ---
    
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    % This operation will automatically promote I to a dlarray
    R = I + sin_theta * K_dl + (1 - cos_theta) * K2;
end