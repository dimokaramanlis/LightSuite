function opts = setupFusiOptions(data, pxsize, opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

opts.sample     = data;
opts.sample_res = pxsize;
%==========================================================================
bofile = fullfile(opts.savepath, 'brain_orientation.txt');
fprintf('Looking for brain orientation data in %s\n', bofile)
if exist(bofile, 'file')
    permvec = load(bofile);
    fprintf('Found it, brain orientation is %s\n', mat2str(permvec))
else
    fprintf('You have to specify the orientation, check GUI\n')
    permvec = getBrainOrientation(data,opts.atlas, [], pxsize);
    % Save to file
    writematrix(permvec, bofile);
    fprintf('Brain orientation saved as %s to %s\n', mat2str(permvec), bofile);
end
%==========================================================================
opts.permute_sample_to_atlas = permvec;
%==========================================================================
save(fullfile(opts.savepath, 'regopts.mat'), '-struct', 'opts')

%==========================================================================
end