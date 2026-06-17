% LS_TRACE_NEUROPIXELS  Trace Neuropixels (or other linear) probe tracks in a
% registered LightSuite brain volume and export AP_histology-style probe_ccf.
%
% Prerequisite: the lightsheet brain must already be registered to the Allen
% CCF (run demos\ls_analyze_lightsheet_volume.m first).  The output folder
% (opts.savepath) must contain regopts.mat and transform_params.mat.

% folder produced by the lightsheet registration pipeline
savepath = 'D:\DATA\DK001\lightsuite';

%% (manual) annotate probe tracks
% Navigate coronal slices, press 1-9 to pick a probe, and click points along
% each probe track (click it in every slice where it appears for a robust
% fit).  Press S to save points, F to fit all probes.
%
% Pressing F transforms the points to atlas space, fits a straight line per
% probe (SVD), reads the brain regions along the trajectory, saves
% <savepath>\probe_ccf.mat, and opens the trajectory figure.
annotateNeuropixelsProbes(savepath);

%% (auto) re-load and re-plot saved results without the GUI
[probe_ccf, cf] = plotNeuropixelsTracingResults(savepath); %#ok<ASGLU>

% probe_ccf is a struct array (one entry per probe) with AP_histology fields:
%   .points            – clicked points in CCF voxels [AP DV ML]
%   .trajectory_coords – [insertion; tip] CCF voxels [AP DV ML]
%   .trajectory_areas  – table of regions along the track (acronym, name,
%                        colour, trajectory_depth [enter exit] µm, n_voxels)
% plus .probe_number and line-fit metadata (.fit_centroid/.fit_direction).

% Example: print the region breakdown for the first probe
disp(probe_ccf(1).trajectory_areas);
