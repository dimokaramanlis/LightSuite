function imgStack = readDownStack(filename, varargin)
%READDOWNSTACK Summary of this function goes here
%   Detailed explanation goes here

% Get information about the file
info = imfinfo(filename);
num_images = numel(info);
num_chans  = numel(info(1).BitsPerSample);

if nargin < 2
    chanread = 1:num_chans;
else
    chanread = varargin{1};
end

% Preallocate the 3D matrix based on the size and type of the first image
imgStack = zeros(info(1).Height, info(1).Width, numel(chanread), ...
    num_images, sprintf('uint%d',info(1).BitsPerSample(1)));

% Read each slice into the 3D matrix
for k = 1:num_images
    currim               = imread(filename, 'Index', k, 'Info', info);
    imgStack(:, :, :, k) = currim(:, :, chanread);
end
imgStack = squeeze(imgStack);

end

