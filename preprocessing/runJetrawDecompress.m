function runJetrawDecompress(pathtarget, pathList)
    % Construct the base command
    
    fprintf('Decompressing tiffs...\n')
    %======================================================================
    Nfiles    = numel(pathList);
    batchsize = min(60, Nfiles);
    Nbatches  = ceil(Nfiles/batchsize);
    %======================================================================
    msg = []; tic;
    for ibatch = 1:Nbatches
        istart = (ibatch - 1)*batchsize + 1;
        iend   = min(ibatch * batchsize, Nfiles);
        commandStr = sprintf('%s %s %s ', 'jetraw decompress -d', pathtarget, pathList{istart:iend});
         % Trim any trailing space
        commandStr = strtrim(commandStr);
        % Run the command
        [status, cmdout] = system(commandStr);
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Batch %d/%d. Time elapsed %2.2f s...\n', ibatch, Nbatches,toc);
        fprintf(msg);
    end
    %======================================================================


   

end