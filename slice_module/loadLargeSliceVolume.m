function imgStack = loadLargeSliceVolume(folderpath, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%-------------------------------------------------------------------------
% Modified code to find both .tif and .tiff files:
files_tif  = dir(fullfile(folderpath, '*.tif'));
files_tiff = dir(fullfile(folderpath, '*.tiff'));
% Concatenate the results from both dir calls
filesget = [files_tif; files_tiff];
filepaths  = fullfile({filesget(:).folder}', {filesget(:).name}');
Nchannels  = numel(filepaths);
%-------------------------------------------------------------------------
if nargin<2
    chanids = 1:Nchannels;
else
    chanids = varargin{1};
end
%-------------------------------------------------------------------------

info       = imfinfo(filepaths{chanids(1)});
num_images = numel(info);

imgStack  = zeros(info(1).Height, info(1).Width, Nchannels, ...
    num_images, sprintf('uint%d',info(1).BitsPerSample(1)));


% Read each slice into the 3D matrix
for ichan = 1:numel(chanids)
    info  = imfinfo(filepaths{chanids(ichan)});
    for k = 1:num_images
        currim                   = imread(filepaths{chanids(ichan)}, 'Index', k, 'Info', info);
        imgStack(:, :, ichan, k) = currim;
    end
end
imgStack = squeeze(imgStack);



end