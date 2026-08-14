%==========================================================================
% get original Allen Atlas
[tvori, avori, parcelinfo, avorires] = loadAtlasInfo('allen2020_10um');

%==========================================================================
%%
atlasread       = load('D:\fusi_test\physio\fUS_3D_images.mat', 'atlas');
permmaceatlas   = [2 1 3];
atlasres        = 50;

tv = atlasread.atlas.Vascular;
tv = permute(tv, permmaceatlas);
av = permute(atlasread.atlas.Regions, permmaceatlas);

tvregori                = imresize3(tvori, 10/atlasres);
avregori                = imresize3(avori, 10/atlasres, "Method","nearest");
%%
dpdata = 'D:\fusi_test\new_atlas';
opts.permute_sample_to_atlas = [1 2 3];
opts.sample = tv;
opts.parcelinfo = parcelinfo;
opts.sample_res = [1 1 1]*0.05;
opts.atlas_res  = [1 1 1]*0.05;
opts.savepath = fullfile(dpdata, 'lightsuite');
opts.mousename = 'atlasfit';
makeNewDir(opts.savepath);
opts.atlas      = tvregori;
opts.annotation = avregori;
%%
matchControlPoints_minimal(opts);


%% (auto) optimize transforms 
opts.bspline_spatial_scale = 2;
multiobjRegistrationFusi(opts, 2, false); %0.2 works


%%
trparams = load(fullfile(opts.savepath, "transform_params.mat"));
Rfix = imref3d(size(tvregori));
Rmov = imref3d(size(tvregori));

paramsfin = elastix_parameter_read(trparams.tform_bspline_samp20um_to_atlas_20um_px);
paramsfin.ResultImagePixelType = 'float';
[savepath, fname, fext] = fileparts(trparams.tform_bspline_samp20um_to_atlas_20um_px);

new_temp_path = fullfile(savepath, sprintf('%s_temp_annotation%s', fname, fext));
elastix_paramStruct2txt(new_temp_path, paramsfin);

% call final transformix
tvnew  = transformix(tv,new_temp_path, ...
    'movingscale', opts.sample_res.*[1 1 1], 'verbose', false);
delete(new_temp_path);

tvnew  = imwarp(tvnew, Rfix, trparams.tform_affine_samp20um_to_atlas_10um_px,'OutputView',Rmov);

tvnew2 = (tvnew(:, :, 1:114) + flip(tvnew(:,:,115:end), 3))/2;
tvnew2 = cat(3, tvnew2, flip(tvnew2,3));
tvnew2 = uint16((2^16-1)*tvnew2/max(tvnew2,[],'all'));

%%
info = struct();
info.ImageSize = size(tvnew2);
info.Datatype           = 'uint16';
info.Description        = '';
info.Version            = '1.0';
info.Qfactor            = 1;
info.PixelDimensions    = opts.atlas_res;
info.SliceCode          = 'Unknown';
info.SpaceUnits         = 'Millimeter';
info.TimeUnits          = 'None';
info.FrequencyDimension = 0;
info.PhaseDimension     = 0;
info.SpatialDimension   = 3;

% save as atlas
niftiwrite(tvnew2,   fullfile('D:\fusi_test\new_atlas\vessel_atlas_template.nii'), info, 'Compressed',true);
% save as annotation
niftiwrite(avregori, fullfile('D:\fusi_test\new_atlas\vessel_atlas_annotation.nii'),  info,'Compressed',true);

%%

% for plotting
movWarpedshow = single(tv);
curlims = quantile(movWarpedshow(movWarpedshow>0), [0.01 0.9],'all');
movWarpedshow = uint8(255 * (movWarpedshow - curlims(1))/range(curlims));
f = figure('Position',[100 100 1000 800]);
ax = gca;
for ii = 1:264; imshowpair(squeeze(tvregori(ii,:,:)), squeeze(movWarpedshow(ii,:,:)),'Parent',ax); title(ii); pause; end

% for ii = 1:264; imagesc(squeeze(bb(ii,:,:))); title(ii); pause; end
