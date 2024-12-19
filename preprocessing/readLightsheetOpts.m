function opts = readLightsheetOpts(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
tim      = dir(fullfile(opts.datafolder, '*.tif'));
tfile    = imread(fullfile(tim(1).folder, tim(1).name));
[Ny, Nx] = size(tfile);
%--------------------------------------------------------------------------
opts.Nz  = numel(tim);
opts.Nx  = Nx; 
opts.Ny  = Ny; 
%--------------------------------------------------------------------------

end