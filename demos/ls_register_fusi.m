% setup atlas to use
[tvvessel, avvessel, parcelinfo, atlas_res] = loadAtlasInfo('allen2020fusi_50um');
%% load anatomy scan

opts.atlas      = tvvessel;
opts.annotation = avvessel;
opts.atlas_res  = atlas_res;
opts.parcelinfo = parcelinfo;
dpdata          = 'D:\fusi_test\physio'; % folder with data
opts.savepath   = fullfile(dpdata, 'lightsuite');
opts.mousename  = 'test_scan';
makeNewDir(opts.savepath);

% load data and pxsize
exampleImage   = load(fullfile(dpdata, 'scan_anatomy.mat'));
pxsizeanatomy  = exampleImage.anatomic.VoxelSize *1e-3;
volumetest     = exampleImage.anatomic.Data;
volumetest     = volumetest/median(volumetest,'all'); % bring values to reasonable range

opts           = setupFusiOptions(volumetest, pxsizeanatomy, opts);
%==========================================================================
%% (manual) add landmarks
matchControlPoints_minimal(opts);
%==========================================================================
%% (auto) optimize transforms 
wtpoints = 0.2;
multiobjRegistrationFusi(opts, wtpoints, false);
%==========================================================================
%% load functional recording
funscan = load('D:\fusi_test\physio\sessionAverage.mat');
scanfus = funscan.scanfus;
recname = 'example_session';
%==========================================================================
%% (auto) bring functional recording to anatomy with rigid transform
tformfuntoanatomy = fusiRecordingToAnatomy(opts, scanfus.Data, scanfus.VoxelSize*1e-3, recname);
%% (auto) apply transform to correlation map
t0 = 30; t1 = 40;
outmap  = mapCorrelation(scanfus, t0, t1);
rescorr = applyFusiTransforms(opts, outmap.Data, outmap.VoxelSize*1e-3, true, tformfuntoanatomy);
save(fullfile(opts.savepath, recname, 'corr_registered.mat'), '-struct', 'rescorr','-v7.3');
%% (auto) apply transform to timeseries
restime = applyFusiTransforms(opts, scanfus.Data, scanfus.VoxelSize*1e-3, false, tformfuntoanatomy);
save(fullfile(opts.savepath, recname, 'time_registered.mat'), '-struct', 'restime','-v7.3');
