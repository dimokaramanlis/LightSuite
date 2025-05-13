function imgStack = loadLargeSliceVolume(folderpath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Modified code to find both .tif and .tiff files:
files_tif  = dir(fullfile(folderpath, '*.tif'));
files_tiff = dir(fullfile(folderpath, '*.tiff'));
% Concatenate the results from both dir calls
filesget = [files_tif; files_tiff];

filepaths  = fullfile({filesget(:).folder}', {filesget(:).name}');
Nchannels  = numel(filepaths);
info       = imfinfo(filepaths{1});
num_images = numel(info);

imgStack  = zeros(info(1).Height, info(1).Width, Nchannels, ...
    num_images, sprintf('uint%d',info(1).BitsPerSample(1)));


% Read each slice into the 3D matrix
for ichan = 1:Nchannels
    info  = imfinfo(filepaths{ichan});
    for k = 1:num_images
        currim                   = imread(filepaths{ichan}, 'Index', k, 'Info', info);
        imgStack(:, :, ichan, k) = currim;
    end
end
imgStack = squeeze(imgStack);



end