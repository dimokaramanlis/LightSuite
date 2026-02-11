% set main data path - where the tiffs are
dpspinesample      = 'D:\spine_registration\sample_test';
bcpdpath           = which('bcpd.exe');
%% load sample and atlas - set resolution and channel for registration
sampleres          = 20; % in micrometers
[cordvol, opts]    = readSpinalCordSample(dpspinesample, sampleres);
%%
opts.regchan       = 2; % choose registration channel
opts.bcpdpath      = bcpdpath;
regopts            = prepareCordSampleForRegistration(cordvol, opts);
%% (manual) trace spinal cord to untwist and unbend
regopts = loadRegOpts(dpspinesample);
spinal_cord_aligner(regopts);

%% (auto) initialize registration
regopts = loadRegOpts(dpspinesample);
initializeCordRegistration(regopts);
%% (manual) add control points
regopts = loadRegOpts(dpspinesample);
matchControlPointsSpine(regopts);

%% (auto) perform nonlinear registration (b-spline)
regopts = loadRegOpts(dpspinesample);
control_point_weight = 0.2;
tparams = multiobjCordRegistration(regopts, control_point_weight);

%% (auto) register all volume channels to the atlas
transform_params = load(fullfile(regopts.lsfolder, 'transform_params.mat'));
regvol = generateRegisteredCordVolume(regopts, transform_params);