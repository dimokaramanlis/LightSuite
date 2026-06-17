function [probe_ccf, cf] = plotNeuropixelsTracingResults(savepath)
% PLOTNEUROPIXELSTRACINGRESULTS  Load and display saved Neuropixels probe
%   tracing results from a LightSuite folder.
%
%   [PROBE_CCF, CF] = plotNeuropixelsTracingResults(SAVEPATH)
%
%   Loads SAVEPATH/probe_ccf.mat (written by annotateNeuropixelsProbes) and
%   re-creates the trajectory figure produced by plotProbeAtlasImages without
%   re-running the annotation GUI.
%
%   Inputs
%   ------
%   SAVEPATH  – LightSuite output folder in which annotateNeuropixelsProbes
%               was run.  Must contain probe_ccf.mat.
%
%   Outputs
%   -------
%   PROBE_CCF – the loaded probe_ccf struct array.
%   CF        – handle to the trajectory figure.

    ccf_file = fullfile(savepath, 'probe_ccf.mat');
    if ~exist(ccf_file, 'file')
        error(['plotNeuropixelsTracingResults: probe_ccf.mat not found in:\n  %s\n' ...
               'Run annotateNeuropixelsProbes and press F to fit probes first.'], savepath);
    end

    s = load(ccf_file, 'probe_ccf');
    if ~isfield(s, 'probe_ccf') || isempty(s.probe_ccf)
        error('plotNeuropixelsTracingResults: %s does not contain a non-empty probe_ccf.', ccf_file);
    end
    probe_ccf = s.probe_ccf;

    fprintf('plotNeuropixelsTracingResults: loaded %d probe(s) from %s\n', ...
        numel(probe_ccf), ccf_file);

    cf = plotProbeAtlasImages(probe_ccf);
    if ~isempty(cf) && isgraphics(cf)
        cf.Name = sprintf('Neuropixels Probe Trajectories – %s', savepath);
    end
end
