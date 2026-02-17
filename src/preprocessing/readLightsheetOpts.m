function opts = readLightsheetOpts(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Looking for data in %s\n', opts.datafolder)
%--------------------------------------------------------------------------
tim      = dir(fullfile(opts.datafolder, '*.tif'));
opts.Nz  = numel(tim);
fprintf('Found %d tiff files ', opts.Nz)
%--------------------------------------------------------------------------
tinfo    = imfinfo(fullfile(tim(1).folder, tim(1).name));
Ny       = tinfo.Height;
Nx       = tinfo.Width;
%--------------------------------------------------------------------------
allnyx = nan(opts.Nz, 2);
for ii = 1:5:opts.Nz
    tinfo         = imfinfo(fullfile(tim(ii).folder, tim(ii).name));
    allnyx(ii, :) = [tinfo.Height tinfo.Width];
end
icheck = ~isnan(sum(allnyx,2));
assert(all(allnyx(icheck, :) == [Ny Nx], "all"), ...
    "some slices do not have matching size, LightSuite cannot proceed")
%--------------------------------------------------------------------------
opts.Nx  = Nx; 
opts.Ny  = Ny; 
fprintf('of size %d x %d px\n', opts.Ny, opts.Nx)
%--------------------------------------------------------------------------

end