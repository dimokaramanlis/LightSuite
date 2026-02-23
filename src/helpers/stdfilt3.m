function [J] = stdfilt3(I, nhoodSize)
% STDFILT3 Robust implementation using separable convolution
%   Fixes "integralImage3" dimension errors caused by NaNs.

    if nargin < 2, nhoodSize = [3 3 3]; end
    
    I = double(I);
    N = prod(nhoodSize); % Total number of elements in the neighborhood
    
    % 1. Pad the image to prevent zero-padding edge artifacts
    % floor(nhoodSize/2) assumes odd-sized dimensions (e.g., 3, 5, 7)
    padSize = floor(nhoodSize / 2);
    I_padded = padarray(I, padSize, 'symmetric');
    
    % 2. Define Separable Kernels
    kH = ones(nhoodSize(1), 1, 1) / nhoodSize(1);
    kW = ones(1, nhoodSize(2), 1) / nhoodSize(2);
    kD = ones(1, 1, nhoodSize(3)) / nhoodSize(3);

    % 3. Helper for Separable Convolution
    % Using 'valid' because we already manually padded the input
    smooth3 = @(V) convn(convn(convn(V, kH, 'valid'), kW, 'valid'), kD, 'valid');

    % 4. Compute Variance: E[X^2] - (E[X])^2
    localMean   = smooth3(I_padded);
    localMeanSq = smooth3(I_padded.^2);
    
    % Ensure no negative variance due to floating point rounding
    localVar = max(0, localMeanSq - localMean.^2);
    
    % 5. Apply Bessel's Correction (N / N-1) for sample variance
    if N > 1
        localVar = localVar * (N / (N - 1));
    end
    
    J = sqrt(localVar);
end