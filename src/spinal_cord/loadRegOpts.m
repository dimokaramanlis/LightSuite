function opts = loadRegOpts(inputpath)
%LOADREGOPTS Load the LightSuite options struct of a sample.
%
%   OPTS = LOADREGOPTS(INPUTPATH) returns the options struct written by the
%   pipeline. INPUTPATH can be
%     * the folder holding regopts.mat (usually OPTS.savepath),
%     * a data folder with a 'lightsuite' subfolder holding it, or
%     * the path to regopts.mat itself.
%
%   Both on-disk layouts are understood: the current one, which saves a single
%   'opts' variable, and the older spinal cord one, which saved the fields of
%   the struct at the top level of the file.
%
%   See also PREPARECORDSAMPLEFORREGISTRATION, INITIALIZECORDREGISTRATION.

if isstruct(inputpath)
    opts = inputpath;
    return
end

inputpath = char(inputpath);

if isfile(inputpath)
    optsfile = inputpath;
else
    candidates = {fullfile(inputpath, 'regopts.mat'), ...
                  fullfile(inputpath, 'lightsuite', 'regopts.mat')};
    ifound     = find(cellfun(@(x) exist(x, 'file') == 2, candidates), 1);
    if isempty(ifound)
        error('loadRegOpts:notFound', ...
            'Could not find regopts.mat in %s (or its lightsuite subfolder).', ...
            inputpath);
    end
    optsfile = candidates{ifound};
end

loaded = load(optsfile);
if isfield(loaded, 'opts')
    opts = loaded.opts;
else
    opts = loaded;   % legacy layout, fields saved at the top level
end

if ~isfield(opts, 'savepath') || isempty(opts.savepath)
    opts.savepath = fileparts(optsfile);
end

end
