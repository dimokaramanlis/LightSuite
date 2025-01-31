function imgStack = readDownStack(filename)
%READDOWNSTACK Summary of this function goes here
%   Detailed explanation goes here

% Get information about the file
info = imfinfo(filename);
num_images = numel(info);

% Preallocate the 3D matrix based on the size and type of the first image
imgStack = zeros(info(1).Height, info(1).Width, num_images, sprintf('uint%d',info(1).BitDepth));

% Read each slice into the 3D matrix
for k = 1:num_images
    imgStack(:, :, k) = imread(filename, 'Index', k, 'Info', info);
end
end

