function saveLightsheetVolume(volsave, savepath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Create a Tiff object
t = Tiff(savepath, 'w');

% Set up TIFF tags
tagstruct.ImageLength = size(volsave, 1);
tagstruct.ImageWidth = size(volsave, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

% Write each slice to the TIFF file
for k = 1:size(volsave, 3)
    setTag(t, tagstruct);
    write(t, volsave(:, :, k));
    if k ~= size(volsave, 3)
        writeDirectory(t);
    end
end

% Close the Tiff object
close(t);

end