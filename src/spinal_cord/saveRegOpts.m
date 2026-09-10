function saveRegOpts(opts)
%SAVEREGOPTS Write the options struct back to regopts.mat in OPTS.savepath.
%
%   SAVEREGOPTS(OPTS) saves OPTS as the single variable 'opts', the layout every
%   LightSuite loader expects. Volumes are kept out of this file - the pipeline
%   stores them as TIFFs and keeps only their paths here - so the plain (v7)
%   format is enough and stays quick to load.
%
%   See also LOADREGOPTS.

assert(isfield(opts, 'savepath') && ~isempty(opts.savepath), ...
    'saveRegOpts:noSavePath', 'opts.savepath must be set before saving.');

makeNewDir(opts.savepath);
save(fullfile(opts.savepath, 'regopts.mat'), 'opts');

end
