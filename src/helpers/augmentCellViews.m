function [dataOut] = augmentCellViews(data)
% AUGMENTCELLVIEWS Applies 3D-consistent augmentation to 3-view cell images.
%
%   Input:
%       data: Cell array {Image, Label} provided by arrayDatastore.
%             Image is [H x W x 3].
%             Channel 1: Y-Z view (Rows=Y, Cols=Z)
%             Channel 2: X-Z view (Rows=X, Cols=Z)
%             Channel 3: X-Y view (Rows=Y, Cols=X)
%
%   Transformations implemented:
%       1. Random Flips along 3D axes (X, Y, Z).
%       2. Random Translations along 3D axes.

    img = data{1};  % Input image [H, W, 3]
    lbl = data{2};  % Label
    
    % --- Settings ---
    maxShift = 0;   % Maximum pixel shift for translation
    probFlip = 0.5; % Probability of flipping an axis
    
    % Get dimensions
    [H, W, ~] = size(img);
    
    % ---------------------------------------------------------
    % 1. 3D REFLECTIONS (Flips)
    % ---------------------------------------------------------
    % We flip the underlying 3D axes, which affects specific image axes per view.
    
    % Flip Dimension 1 (Y-axis): Affects Rows of Ch1 and Rows of Ch3
    if rand() < probFlip
        img(:,:,1) = flipud(img(:,:,1)); 
        img(:,:,3) = flipud(img(:,:,3));
    end
    
    % Flip Dimension 2 (X-axis): Affects Rows of Ch2 and Cols of Ch3
    if rand() < probFlip
        img(:,:,2) = flipud(img(:,:,2));
        img(:,:,3) = fliplr(img(:,:,3));
    end
    
    % Flip Dimension 3 (Z-axis): Affects Cols of Ch1 and Cols of Ch2
    if rand() < probFlip
        img(:,:,1) = fliplr(img(:,:,1));
        img(:,:,2) = fliplr(img(:,:,2));
    end
    
    % ---------------------------------------------------------
    % 2. 3D TRANSLATIONS (Shifts)
    % ---------------------------------------------------------
    % Generate random integer shifts for Y, X, Z
    % shifts = randi([-maxShift, maxShift], 1, 3); 
    % dy = shifts(1); % Shift Y
    % dx = shifts(2); % Shift X
    % dz = shifts(3); % Shift Z
    % 
    % % Apply coupled shifts:
    % % View 1 (Y-Z): Shift Y (Rows), Shift Z (Cols)
    % img(:,:,1) = imtranslate(img(:,:,1), [dz, dy], 'FillValues', 0);
    % 
    % % View 2 (X-Z): Shift X (Rows), Shift Z (Cols)
    % % Note: In Ch2, X is on Rows. A "vertical" image shift corresponds to X.
    % img(:,:,2) = imtranslate(img(:,:,2), [dz, dx], 'FillValues', 0);
    % 
    % % View 3 (X-Y): Shift X (Cols), Shift Y (Rows)
    % img(:,:,3) = imtranslate(img(:,:,3), [dx, dy], 'FillValues', 0);

    dataOut = {img, lbl};
end