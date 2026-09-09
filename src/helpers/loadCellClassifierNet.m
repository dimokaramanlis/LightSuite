function [net, netinfo] = loadCellClassifierNet(netin)
%LOADCELLCLASSIFIERNET Resolve a cell-classifier network from a path or object.
%   [NET, NETINFO] = LOADCELLCLASSIFIERNET(NETIN) accepts any of the forms the
%   pipeline may be handed a trained artifact classifier in and returns the
%   network object itself:
%
%     * a char/string path to a .mat file saved by demos/trainClassificationNetwork.m
%       (a 'net' variable, optionally alongside the accuracies);
%     * a struct holding a 'net' field (e.g. the output of load());
%     * the network object itself (SeriesNetwork, DAGNetwork or dlnetwork).
%
%   NETINFO records where the network came from, so the classification results
%   written to disk can say which network produced them:
%     .source     'file' | 'struct' | 'object'
%     .path       the .mat path, or '' when a network object was passed in
%     .name       a short label for messages
%     .accuracy   validation accuracy if the file carried one, else NaN
%
%   See also CLASSIFYDETECTEDCELLS, PREPAREIMAGESFORCNN.

netinfo = struct('source', 'object', 'path', '', 'name', '', 'accuracy', NaN);

if ischar(netin) || isstring(netin)
    netpath = char(netin);
    if ~isfile(netpath)
        error('loadCellClassifierNet:fileNotFound', ...
            'Classifier network file not found: %s', netpath);
    end
    dat = load(netpath);
    [~, nname] = fileparts(netpath);
    netinfo.source = 'file';
    netinfo.path   = netpath;
    netinfo.name   = nname;
    if isfield(dat, 'valaccuracy')
        netinfo.accuracy = dat.valaccuracy;
    end
    net = netFromStruct(dat, netpath);

elseif isstruct(netin)
    netinfo.source = 'struct';
    netinfo.name   = 'struct';
    if isfield(netin, 'valaccuracy')
        netinfo.accuracy = netin.valaccuracy;
    end
    net = netFromStruct(netin, 'the supplied struct');

else
    net          = netin;
    netinfo.name = class(netin);
end

if ~(isa(net, 'SeriesNetwork') || isa(net, 'DAGNetwork') || isa(net, 'dlnetwork'))
    error('loadCellClassifierNet:badNetwork', ...
        ['Expected a SeriesNetwork, DAGNetwork or dlnetwork, got %s. Pass the ' ...
         'network trained by demos/trainClassificationNetwork.m, or the path ' ...
         'of the .mat file it was saved to.'], class(net));
end

end

%==========================================================================
function net = netFromStruct(dat, srcname)
%NETFROMSTRUCT Pull the network out of a loaded .mat struct.
if isfield(dat, 'net')
    net = dat.net;
    return;
end

% no 'net' field: accept the file if exactly one variable looks like a network
fn     = fieldnames(dat);
isnet  = false(numel(fn), 1);
for k = 1:numel(fn)
    v = dat.(fn{k});
    isnet(k) = isa(v, 'SeriesNetwork') || isa(v, 'DAGNetwork') || isa(v, 'dlnetwork');
end

if nnz(isnet) == 1
    net = dat.(fn{isnet});
elseif nnz(isnet) == 0
    error('loadCellClassifierNet:noNet', ...
        ['%s holds no network. Expected a ''net'' variable, as saved by ' ...
         'demos/trainClassificationNetwork.m.'], srcname);
else
    error('loadCellClassifierNet:ambiguousNet', ...
        ['%s holds %d networks (%s) and no ''net'' variable, so it is ' ...
         'ambiguous which one to use.'], srcname, nnz(isnet), ...
        strjoin(fn(isnet), ', '));
end
end
