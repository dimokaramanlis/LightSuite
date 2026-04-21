function visualizeCellDetections(inputpath,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) isstring(x) || ischar(x));
validFinishes = {'sample','atlas'};
checkFinish   = @(x) any(validatestring(x,validFinishes));
addParameter(p, 'Space', 'sample',checkFinish);
addParameter(p, 'Maxpoints', 8e4, @isnumeric);
parse(p, inputpath, varargin{:});
params = p.Results;
%==========================================================================
issamp   = false;
opts     = load(fullfile(inputpath, "regopts.mat"));
switch params.Space
    case 'sample'
        lastbit = "*_locations_sample.mat";
        issamp  = true;
        pxsize = opts.opts.pxsize;
        txtuse = 'Sample space';
    case 'atlas'
        lastbit = "*_locations_atlas.mat";
        pxsize = opts.opts.atlasres*[1 1 1];
        txtuse = 'Atlas space';
end

locpaths = dir(fullfile(inputpath, lastbit));
%==========================================================================
if isempty(locpaths)
    warning("Can't find any valid detections...");
    return
end
%==========================================================================
Nchans = numel(locpaths);
fh = figure('Position',[50, 50, 500 + 500*(Nchans-1), 700]);
p  = panel();
p.pack('h', numel(locpaths));
for ipath = 1:numel(locpaths)
    dpcurr = fullfile(locpaths(ipath).folder, locpaths(ipath).name);
    if issamp
        clocs  = load(dpcurr, "cell_locations");
        clocs  = clocs.cell_locations;
    else
        clocs  = load(dpcurr, "atlasptcoords");
        clocs  = clocs.atlasptcoords;
    end
    clocs = subsampplot(clocs(:, 1:3), params.Maxpoints);
    
    axcurr = p(ipath).select();
    if issamp
        clocs = clocs.*pxsize;
        scatter3(clocs(:,1), clocs(:,2), clocs(:,3), 2, 'filled');
    else
        plotBrainGrid([], axcurr); hold on;
        scatter3(clocs(:,2), clocs(:,3), clocs(:,1), 2, 'filled');
    end
    axis equal; axis tight;
    ax = gca; ax.ZDir = 'reverse';
    nameuse = strrep(locpaths(ipath).name,'_',' ');

    title(sprintf('%s, %s', txtuse, nameuse))
end

%==========================================================================
end

% plotBrainGrid; hold on;
% scatter3(atlasptcoords(iplot,2),atlasptcoords(iplot,3),atlasptcoords(iplot,1),2,...
%     'filled','MarkerFaceAlpha',0.5)