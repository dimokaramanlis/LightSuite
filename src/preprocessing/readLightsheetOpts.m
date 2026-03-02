function opts = readLightsheetOpts(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Looking for data in %s\n', opts.datafolder)
%--------------------------------------------------------------------------
tifftype = getOr(opts, 'tifftype', 'planeperfile');
tiffiles = dir(fullfile(opts.datafolder, '*.tiff'));
tifiles  = dir(fullfile(opts.datafolder, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
opts.multitiffs = false;
%--------------------------------------------------------------------------
switch tifftype
    %======================================================================
    case 'planeperfile' % classic for Terastitcher
        fprintf('Using planeperfile loading\n')
        %--------------------------------------------------------------------------
        opts.Nz  = numel(tifiles);
        fprintf('Found %d tiff files ', opts.Nz)
        %--------------------------------------------------------------------------
        tinfo    = imfinfo(fullfile(tifiles(1).folder, tifiles(1).name));
        Ny       = tinfo.Height;
        Nx       = tinfo.Width;
        %--------------------------------------------------------------------------
        allnyx = nan(opts.Nz, 2);
        for ii = 1:5:opts.Nz
            tinfo         = imfinfo(fullfile(tifiles(ii).folder, tifiles(ii).name));
            allnyx(ii, :) = [tinfo.Height tinfo.Width];
        end
        icheck = ~isnan(sum(allnyx,2));
        assert(all(allnyx(icheck, :) == [Ny Nx], "all"), ...
            "some slices do not have matching size, LightSuite cannot proceed")
        %--------------------------------------------------------------------------
        opts.Nx     = Nx; 
        opts.Ny     = Ny; 
        opts.Nchans = 1;
        fprintf('of size %d x %d px\n', opts.Ny, opts.Nx)
        %======================================================================
    case 'channelperfile' % classic for BigStitcher
        fprintf('Using channelperfile loading, every channel should be in a different tiff\n')
        Nfiles = numel(tfiles);
        if Nfiles == 1
             datainfo = BioformatsImage(fullfile(tfiles.folder, tfiles.name));
             opts.Nchans = datainfo.sizeC;
             fprintf('Found a single tiff with %d channels \n', opts.Nchans);
             allnyxz = [datainfo.height datainfo.width datainfo.sizeZ];
        else
            opts.Nchans  = numel(tfiles);
            fprintf('Assuming each channel is a separate tiff. Found %d channels \n', opts.Nchans)
            allnyxz = nan(opts.Nchans, 3);
            for ichan = 1 : opts.Nchans
                datainfo          = BioformatsImage(fullfile(tfiles(ichan).folder, tfiles(ichan).name));
                allnyxz(ichan, :) = [datainfo.height datainfo.width datainfo.sizeZ];
            end
            opts.multitiffs = true;
            assert(all(allnyxz == allnyxz(1,:), "all"), ...
                "some volumes do not have matching size, LightSuite cannot proceed")
        end
        %--------------------------------------------------------------------------
        opts.Nx  = allnyxz(1, 2); 
        opts.Ny  = allnyxz(1, 1); 
        opts.Nz  = allnyxz(1, 3); 
        fprintf('Volume size is %d x %d x %d px\n', opts.Ny, opts.Nx,  opts.Nz)
        %======================================================================
end
%--------------------------------------------------------------------------
opts.tfiles = fullfile({tfiles(:).folder},{tfiles(:).name})';
makeNewDir(opts.savepath)
%--------------------------------------------------------------------------

end