function res = classifyDetectedCells(net, cell_images, imwindow, opts)
%CLASSIFYDETECTEDCELLS Run the CNN artifact classifier over detected candidates.
%   RES = CLASSIFYDETECTEDCELLS(NET, CELL_IMAGES, IMWINDOW) applies the trained
%   cell/noise classifier (demos/trainClassificationNetwork.m) to every detected
%   candidate and returns which of them are real cells. CELL_IMAGES is the
%   [Ncells x Nfeatures] array of flattened three-view cell images written by
%   the detector when opts.savecellimages is on, and IMWINDOW is its half-window
%   (saved next to it; 6*[3 3 2] for the default detector settings).
%
%   The classifier sees the same three maximum-intensity views used for training,
%   built by prepareImagesForCNN.
%
%   MEMORY. A full brain can hold more than a million candidates, and the
%   prepared images are ~16 kB each, so they are never all materialized at once:
%   cells are processed in chunks of opts.batchsize. prepareImagesForCNN
%   rescales its output by a min/max taken over everything it is handed, which
%   would make a chunked result depend on the chunk boundaries, so the constants
%   are measured once over the whole file in a first pass and then reused for
%   every chunk. That makes the chunked result identical to classifying the file
%   in one go.
%
%   Input:
%     net          trained network, or anything loadCellClassifierNet accepts.
%     cell_images  [Ncells x Nfeatures] flattened cell views.
%     imwindow     1x3 half-window of the views.
%     opts         (optional) struct:
%        .batchsize      candidates prepared at once (default 10000).
%        .minibatchsize  network mini-batch (default 256).
%        .goodclass      class counted as a real cell. Default: the class named
%                        '1' (CellLabelingTool writes 1 = Cell, 0 = Noise).
%        .verbose        default true.
%
%   Output struct res:
%     .isgood      [Ncells x 1] logical, true for candidates kept as cells
%     .labels      [Ncells x 1] categorical network output
%     .scores      [Ncells x Nclasses] single class scores
%     .classnames  the network's classes, in score-column order
%     .goodclass   the class treated as "cell"
%     .ncells      Ncells
%     .fracgood    mean(isgood)
%
%   See also PREPAREIMAGESFORCNN, LOADCELLCLASSIFIERNET, TRANSFORMPOINTSTOATLAS.

if nargin < 4, opts = struct(); end
batchsize     = getOr(opts, 'batchsize',     10000);
minibatchsize = getOr(opts, 'minibatchsize', 256);
goodclass     = getOr(opts, 'goodclass',     []);
verbose       = getOr(opts, 'verbose',       true);

net = loadCellClassifierNet(net);

if isempty(cell_images)
    error('classifyDetectedCells:noImages', ...
        ['No cell images to classify. The classifier needs the ''cell_images'' ' ...
         'array, which the detector only writes when opts.savecellimages is true.']);
end

imwindow = double(imwindow(:)).';
if numel(imwindow) ~= 3
    error('classifyDetectedCells:badWindow', 'imwindow must be a 1x3 vector.');
end

Ncells     = size(cell_images, 1);
classnames = classesOf(net);
goodclass  = resolveGoodClass(classnames, goodclass);

edges  = [1:batchsize:Ncells, Ncells + 1];
nbatch = numel(edges) - 1;

%--------------------------------------------------------------------------
% pass 1: normalization constants over the whole file (see note above)
%--------------------------------------------------------------------------
nrm = [inf, -inf];
if verbose
    fprintf('    measuring image scale over %d candidates...\n', Ncells);
end
for ib = 1:nbatch
    idx = edges(ib):edges(ib+1) - 1;
    [~, thisnrm] = prepareImagesForCNN(cell_images(idx, :), imwindow);
    nrm = [min(nrm(1), thisnrm(1)), max(nrm(2), thisnrm(2))];
end

%--------------------------------------------------------------------------
% pass 2: classify
%--------------------------------------------------------------------------
labels = categorical(repmat("", Ncells, 1), cellstr(classnames));
scores = zeros(Ncells, numel(classnames), 'single');

msg = ''; tic;
for ib = 1:nbatch
    idx = edges(ib):edges(ib+1) - 1;
    X   = prepareImagesForCNN(cell_images(idx, :), imwindow, nrm);

    if isa(net, 'dlnetwork')
        sc = predict(net, dlarray(single(X), 'SSCB'));
        sc = single(extractdata(sc)).';
        [~, imax]   = max(sc, [], 2);
        labels(idx) = classnames(imax);
    else
        [lb, sc]    = classify(net, X, 'MiniBatchSize', minibatchsize);
        labels(idx) = lb;
    end
    scores(idx, :) = single(sc);

    if verbose
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('    classified %d/%d candidates. Time elapsed %2.2f s\n', ...
            idx(end), Ncells, toc);
        fprintf(msg);
    end
end

res            = struct();
res.isgood     = labels == goodclass;
res.labels     = labels;
res.scores     = scores;
res.classnames = classnames;
res.goodclass  = goodclass;
res.ncells     = Ncells;
res.fracgood   = mean(res.isgood);

if verbose
    fprintf('    kept %d/%d candidates as cells (%2.1f%%).\n', ...
        nnz(res.isgood), Ncells, 100*res.fracgood);
end

end

%==========================================================================
function cn = classesOf(net)
%CLASSESOF The network's output classes, as a categorical row.
cn = [];
if isprop(net, 'Layers') && ~isempty(net.Layers)
    last = net.Layers(end);
    if isprop(last, 'Classes') && ~isempty(last.Classes)
        cn = last.Classes(:).';
    end
end
if isempty(cn)
    error('classifyDetectedCells:noClasses', ...
        ['Cannot read the class names off this network. Retrain with a ' ...
         'classificationLayer, or classify the cells yourself and pass the ' ...
         'result in.']);
end
end

%--------------------------------------------------------------------------
function gc = resolveGoodClass(classnames, requested)
%RESOLVEGOODCLASS Which class means "this is a real cell".
if ~isempty(requested)
    gc = categorical(string(requested), cellstr(classnames));
    if ismissing(gc)
        error('classifyDetectedCells:badGoodClass', ...
            '''goodclass'' %s is not one of the network''s classes (%s).', ...
            string(requested), strjoin(cellstr(classnames), ', '));
    end
    return;
end

names = string(classnames);

% CellLabelingTool labels 1 = Cell, 0 = Noise, so a "1" class is the cells
ihit = find(names == "1", 1);
if isempty(ihit)
    ihit = find(strcmpi(names, "cell"), 1);
end
if isempty(ihit)
    error('classifyDetectedCells:ambiguousGoodClass', ...
        ['Cannot tell which of the network''s classes (%s) means "real cell". ' ...
         'Pass it explicitly, e.g. ''goodclass'', "%s".'], ...
        strjoin(cellstr(classnames), ', '), names(end));
end
gc = classnames(ihit);
end
