function [meanMap, varMap] = atlasCellMapVariance(mousecells, avsize, sigma, varargin)
% atlasCellMapVariance Calculates weighted mean and variance of cell densities.
%
%   Returns:
%   meanMap: The weighted average density map (same as your original output)
%   varMap:  The weighted variance map across mice

    % --- 1. Setup Inputs ---
    if iscell(mousecells)
        Nmice = numel(mousecells);
    else
        Nmice = 1;
        mousecells = {mousecells};
    end

    % --- 2. Setup Weights ---
    if nargin < 4
        wts = ones(Nmice, 1) / Nmice;
    else
        wts = varargin{1}(:);
        % Normalize weights to sum to 1 immediately to simplify math
        wts = wts / sum(wts);
    end

    % --- 3. Initialize Accumulators (Single Precision for Memory) ---
    % M1 = Weighted Sum (First Moment) -> E[X]
    % M2 = Weighted Sum of Squares (Second Moment) -> E[X^2]
    Nmid = avsize(3)/2;
    avsizeset = avsize;
    avsizeset(3) = Nmid;
    M1 = zeros(avsizeset, 'single');
    M2 = zeros(avsizeset, 'single');
    

    % --- 4. The "One Pass" Loop ---
    for ii = 1:Nmice
        % A. Build Sparse Grid for THIS mouse
        % Note: We use 1.0 as the value here, we apply weights later
        % to the full smoothed volume to ensure correct density scaling.
        thisMouseSparse = ndSparse.build(double(mousecells{ii}(:,[2 1 3])), 1, avsize);
        
        % B. Convert to Full Single
        thisVol = single(full(thisMouseSparse));
        
        % C. Symmetrize (Must happen per mouse now)
        thisVol = (thisVol(:,:, 1:Nmid) + flip(thisVol(:,:, (Nmid + 1):end), 3)) / 2;
        
        % D. Gaussian Smoothing (Must happen before squaring)
        % This represents the "Density Map" (x_i) for this mouse
        smoothedVol = imgaussfilt3(thisVol, sigma);
        
        % E. Accumulate Moments
        % M1 += w_i * x_i
        M1 = M1 + (wts(ii) * smoothedVol);
        
        % M2 += w_i * (x_i^2)
        M2 = M2 + (wts(ii) * (smoothedVol .^ 2));
    end

    % --- 5. Final Calculation ---
    
    % The Weighted Mean (Same as your original mousemap)
    meanMap = M1; 
    
    % The Biased Weighted Variance: E[X^2] - (E[X])^2
    varMap_biased = M2 - (M1 .^ 2);
    
    % Numerical stability check: Variance cannot be negative, but floating point 
    % errors can cause tiny negative numbers (e.g., -1e-10). Clamp to 0.
    varMap_biased(varMap_biased < 0) = 0;

    % --- 6. Unbiased Correction (Optional but Recommended) ---
    % Since your weights sum to 1, we use the reliability weight correction:
    % Factor = 1 / (1 - sum(w^2))
    correctionFactor = 1 / (1 - sum(wts.^2));
    
    varMap = varMap_biased * correctionFactor;
    
end