function [keptpts, classres] = applyCellClassifierCached(inputpts, ds, net, netinfo, ...
    savepath, reclassify)
%APPLYCELLCLASSIFIERCACHED Keep only the candidates the CNN accepts as real cells.
%
%   [KEPTPTS, CLASSRES] = APPLYCELLCLASSIFIERCACHED(INPUTPTS, DS, NET, NETINFO,
%   SAVEPATH, RECLASSIFY) runs the trained artifact classifier NET on the cell
%   images carried by the dataset DS (see READCELLLOCATIONSFILE) and returns the
%   subset of INPUTPTS the network accepts.
%
%   The result is cached as '<ds.outbase>_classification.mat' in SAVEPATH, so
%   the network runs once per point file and later runs reuse it. Pass
%   RECLASSIFY true to run it again anyway; a cache that no longer covers the
%   same number of candidates is discarded automatically.
%
%   Point sets that carry no cell images - CSV and XML files, and .mat
%   detections saved with opts.savecellimages off - are returned unfiltered with
%   a warning, since there is nothing for the network to look at.
%
%   Shared by TRANSFORMPOINTSTOATLAS and TRANSFORMCORDPOINTSTOATLAS.
%
%   See also CLASSIFYDETECTEDCELLS, LOADCELLCLASSIFIERNET, READCELLLOCATIONSFILE.

classres = [];
keptpts  = inputpts;

if isempty(ds.cell_images)
    warning('applyCellClassifierCached:noCellImages', ...
        ['%s carries no cell images, so the classifier cannot be run on it; ' ...
         'its points are transformed unfiltered. Cell images are only written ' ...
         'for .mat detections saved with opts.savecellimages = true.'], ds.label);
    return;
end

cachefile = fullfile(savepath, sprintf('%s_classification.mat', ds.outbase));
ncells    = size(inputpts, 1);

if exist(cachefile, 'file') && ~reclassify
    cached = load(cachefile);
    if isfield(cached, 'isgood') && numel(cached.isgood) == ncells
        fprintf('  Reusing cached classification (%s).\n', cachefile);
        classres = cached;
        keptpts  = inputpts(cached.isgood, :);
        return;
    end
    warning('applyCellClassifierCached:staleClassification', ...
        ['Cached classification %s covers %d candidates but %s now has %d; ' ...
         're-running the network.'], cachefile, ...
        numel(getOr(cached, 'isgood', [])), ds.label, ncells);
end

fprintf('  Classifying %d candidates with %s...\n', ncells, netinfo.name);
res = classifyDetectedCells(net, ds.cell_images, ds.imwindow);

classres             = res;
classres.netname     = netinfo.name;
classres.netpath     = netinfo.path;
classres.netaccuracy = netinfo.accuracy;
classres.sourcefile  = ds.sourcefile;
classres.created     = datetime('now');

save(cachefile, '-struct', 'classres');
fprintf('  Saved classification to %s\n', cachefile);

keptpts = inputpts(res.isgood, :);

end
