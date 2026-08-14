function res = fusiRecordingToAnatomy(opts, data, voxelsize_mm, recname, showimage)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    recname = '';
end

if nargin < 5
    showimage = true;
end
outpath    = fullfile(opts.savepath, recname);
makeNewDir(outpath);

permvec =  opts.permute_sample_to_atlas;
%---------------------------------------------------------------------------
% prepare movie
moviereg = median(data, 4);
moviereg = moviereg./median(data,'all'); % bring to reasonable scale
moviereg = permuteBrainVolume(moviereg, permvec);
facmovie = voxelsize_mm(abs(permvec))./opts.atlas_res;
moviereg = imresize3(moviereg, 'Scale', facmovie);
Rmovie   = imref3d(size(moviereg));
%---------------------------------------------------------------------------
% prepare anatomy
volumereg  = permuteBrainVolume(opts.sample, permvec);
volumereg  = imresize3(volumereg, 'Scale', opts.sample_res(abs(permvec))./opts.atlas_res);
Ranatomy   = imref3d(size(volumereg));
%---------------------------------------------------------------------------
targetfile = fullfile(outpath,"rigid_atlas_to_samp_20um.txt");
if exist(targetfile, "file")
    fprintf('Found previous rigid transform, loading...\n')
else
    fprintf('Fitting rigid transform.../n')
    [~,~, targetfile, ~] = performElastixRigidRegistration(volumereg, moviereg, 1,outpath);
end
tformrigid              = parse_elastix_tform_all(targetfile);
tformrigid              = rigidtform3d(tformrigid.A);
tformrigid              = tformrigid.invert();

res = struct();
res.tformrigid = tformrigid;
res.Ranatomy   = Ranatomy;
%---------------------------------------------------------------------------
% get some diagnostics to check
outvolel      = imwarp(moviereg, Rmovie, tformrigid, 'OutputView',Ranatomy);
moviequants   = quantile(moviereg, [0.01 0.95],'all');
movietoshow   = uint8(255*(moviereg-moviequants(1))/range(moviequants));
outvolshow    = uint8(255*(outvolel-moviequants(1))/range(moviequants));
anatomyquants = quantile(volumereg, [0.01 0.95],'all');
anatomyshow   = uint8(255*(volumereg-anatomyquants(1))/range(anatomyquants));

if showimage
    cf = visualizeTransformMatch(anatomyshow, cat(4, movietoshow, outvolshow), 1);
    print(cf, fullfile(outpath, sprintf('dim%d_rigid_registration', 1)), '-dpng')
    close(cf);
end

res.anatomyshow = anatomyshow;
res.movietoshow = movietoshow;
res.outvolshow  = outvolshow;
%---------------------------------------------------------------------------
end