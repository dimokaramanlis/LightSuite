function finalvol = getRigidlyAlignedVolume(sliceinfo, slicevol, tformslices, atlasframe);
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[Nchan, Nslices] = size(slicevol, [3 4]);
pxsamp        = sliceinfo.px_process/sliceinfo.px_register;
pxatlas       = sliceinfo.px_atlas/sliceinfo.px_register;
% R_out_fullres = imref2d(atlasframe, pxatlas, pxatlas);

R_out_fullres = imref2d(ceil(atlasframe*pxatlas/pxsamp), pxsamp, pxsamp);

Rsample       = imref2d(size(slicevol,[1 2]), pxsamp, pxsamp);
finalvol      = zeros([R_out_fullres.ImageSize Nchan Nslices], 'uint16');
slicetic = tic; msg = [];
for islice = 1:Nslices

    % we find the x-y points in the original atlas space

    tform = tformslices(islice);
    % tform.Translation = tform.Translation;
    % fval  = sliceinfo.backvalues(1, islice);
    % testslice =  imwarp(slicevol(:, :, 1, islice), Rsample, tform, ...
    %             'linear', 'OutputView', R_out_fullres, 'FillValues', fval);
    for ichan = 1:Nchan
        fval  = sliceinfo.backvalues(ichan, islice);
        finalvol(:,:,ichan, islice) =  imwarp(slicevol(:, :, ichan, islice), Rsample, tform, ...
                'linear', 'OutputView', R_out_fullres, 'FillValues',fval);
    end
     fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d done. Took %2.2f s\n', islice, Nslices, toc(slicetic));
    fprintf(msg);

end

end