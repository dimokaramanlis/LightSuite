function [J] = stdfilt3(I, nhoodSize)
% STDFILT3 Robust implementation using separable convolution
%   Fixes "integralImage3" dimension errors caused by NaNs.

    if nargin < 2, nhoodSize = [3 3 3]; end
    
    % 1. Handle NaNs and Convert to Double
    I = double(I);
    
    % Optional: Fill NaNs with 0 (or nearest) to prevent NaN propagation
    % If you want NaNs to propagate (output = NaN where input = NaN), skip this.
    % maskNaN = isnan(I);
    % I(maskNaN) = 0; 

    % 2. Define Separable Kernels (Box Filter)
    % This is mathematically identical to imboxfilt3 but uses convn
    kH = ones(nhoodSize(1), 1, 1) / nhoodSize(1);
    kW = ones(1, nhoodSize(2), 1) / nhoodSize(2);
    kD = ones(1, 1, nhoodSize(3)) / nhoodSize(3);

    % 3. Helper for Separable Convolution
    % Convolves in Y, then X, then Z sequentially. 
    % 'same' padding keeps output size equal to input size.
    smooth3 = @(V) convn(convn(convn(V, kH, 'same'), kW, 'same'), kD, 'same');

    % 4. Compute Variance: E[X^2] - (E[X])^2
    localMean   = smooth3(I);
    localMeanSq = smooth3(I.^2);
    
    % Ensure no negative variance due to floating point rounding
    localVar = max(0, localMeanSq - localMean.^2);
    
    J = sqrt(localVar);
end