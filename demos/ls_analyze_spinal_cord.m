% set main data path - where the tiffs are
dpspinesample      = 'D:\spine_registration\sample1';
bcpdpath           = 'C:\Users\karamanl\Documents\GitHub\bcpd\win\bcpd.exe';
%% load sample and atlas - set resolution and channel for registration
sampleres          = 20; % in micrometers
[cordvol, opts]    = readSpinalCordSample(dpspinesample, sampleres);
%%
opts.dpspineatlas  = 'D:\AllenAtlas\extra_spine\allen_cord_20um_v1.1';
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
tparams = multiobjCordRegistration(regopts, 0.2);

%% (auto) register all volume channels to the atlas
transform_params = load(fullfile(regopts.lsfolder, 'transform_params.mat'));
regvol = generateRegisteredCordVolume(regopts, transform_params);