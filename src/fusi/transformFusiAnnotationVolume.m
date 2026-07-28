function [volumereg, avnonrigid, avaffine] = transformFusiAnnotationVolume(opts)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

transformparams = fullfile(opts.savepath,"transform_params.mat");
trstruct        = load(transformparams);
permvec         = opts.permute_sample_to_atlas;


volumereg  = permuteBrainVolume(opts.sample, permvec);
volumereg  = imresize3(volumereg, 'Scale', opts.sample_res(abs(permvec))./opts.atlas_res);
Ranatomy   = imref3d(size(volumereg));
Ratlas     = imref3d(size(opts.atlas));
av         = resampleAnnotation(opts.annotation, opts.parcelinfo, 'structure');

avaffine = imwarp(av, Ratlas, trstruct.tform_affine_samp20um_to_atlas_10um_px.invert,...
    'nearest','OutputView',Ranatomy);
bspltformpath = fullfile(opts.savepath,'bspline_atlas_to_samp_20um.txt');
avnonrigid = transformAnnotationVolume(bspltformpath, avaffine, opts.atlas_res);

end