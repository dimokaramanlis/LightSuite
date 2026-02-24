function opts = preprocessLightSheetVolume(opts)
%PREPROCESSLIGHTSHEETVOLUME Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% other info
scaledownxy  = opts.pxsize(1)/opts.registres;
scaledownz   = opts.pxsize(3)/opts.registres;
[Ny, Nx, Nz] = deal(opts.Ny, opts.Nx,opts.Nz);
%--------------------------------------------------------------------------
% options for tiff saving
saveopts.compress  = 'lzw';
saveopts.message   = false;
%--------------------------------------------------------------------------
regvolpaths     = cell(opts.Nchans, 1);
%--------------------------------------------------------------------------
for ichannel = 1:opts.Nchans
    %----------------------------------------------------------------------
    hascells = ismember(ichannel, opts.channelforcells);
    %----------------------------------------------------------------------
    % initialize collections
    backvol  = zeros(ceil(scaledownxy*Ny), ceil(scaledownxy*Nx), Nz, 'uint16');
    if hascells
        fproc    = fullfile(opts.fproc, sprintf('chan_%d_binary_%s.dat', ichannel, opts.mousename));
        fid      = fopen(fproc, 'W');
    end
    msg = []; proctic = tic;
    %----------------------------------------------------------------------
    if ~contains(opts.tifftype, 'planeperfile')
        currbf = BioformatsImage(opts.tfiles{ichannel});
    end
    %----------------------------------------------------------------------
    for islice = 1:Nz
        switch opts.tifftype
            case 'channelperfile'
                currim = currbf.getPlane(islice,1,1);
            case 'planeperfile'
                currim = imread(opts.tfiles{islice});
        end

        % write data
        backvol(:, :, islice) = imresize(currim, scaledownxy);
        
        if hascells
            % we write a median-filtered version to remove artifacts
            fwrite(fid, medfilt2(currim, [3 3], 'symmetric'), "uint16");
        end
        %----------------------------------------------------------------------
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Channel %d/%d. Slice %d/%d. Time per slice %2.2f s. Time elapsed %2.2f s...\n',...
            ichannel, opts.Nchans, islice, Nz, toc(proctic)/islice, toc(proctic));
        fprintf(msg);
        %----------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    fprintf('Saving volume for registration... ')
    voldown             = imresize3(backvol, 'Scale', [1 1 scaledownz]);
    % save volume for control point and registration
    samplepath = fullfile(opts.savepath, sprintf('chan_%d_sample_register_%dum.tif', ...
        ichannel, opts.registres));
  
    if exist(samplepath, 'file')
        delete(samplepath);
    end
    saveastiff(voldown, samplepath, saveopts);
    fprintf('Done!\n')
    %----------------------------------------------------------------------
    if hascells
        %----------------------------------------------------------------------
        fclose(fproc);
        %----------------------------------------------------------------------
        fprintf('Extracting cell candidates from channel %d \n', ichannel)
        opts.prefix     = sprintf('chan_%d_', ichannel);
        peakvalsextract = extractCellsFromVolumeNew(fproc, opts);
        %----------------------------------------------------------------------
        delete(opts.fproc);
        %----------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
end
opts.regvolpath = regvolpaths{opts.channelforregister};
save(fullfile(opts.savepath, sprintf('regopts.mat')), 'opts')
%--------------------------------------------------------------------------

end