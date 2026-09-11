function visualizeCellDetections(inputpath,varargin)
%VISUALIZECELLDETECTIONS Plot the cell detections of every channel in 3D.
%
%   VISUALIZECELLDETECTIONS(SAVEPATH) plots the detections in the LightSuite
%   folder SAVEPATH in sample space, one panel per channel. SAVEPATH can also
%   be the procpath of a slice volume.
%
%   VISUALIZECELLDETECTIONS(SAVEPATH, 'Space', 'atlas') plots the atlas-space
%   detections instead, on the outline of the atlas: PLOTBRAINGRID for points
%   from TRANSFORMPOINTSTOATLAS, PLOTSPINALCORDGRID for points from
%   TRANSFORMCORDPOINTSTOATLAS. Which one applies is read from the files. Cord
%   points the transform extrapolated past the atlas volume are left out.
%
%   'Maxpoints' (default 8e4) caps the number of points drawn per channel.
%
%   See also PLOTBRAINGRID, PLOTSPINALCORDGRID.
%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) isstring(x) || ischar(x));
validFinishes = {'sample','atlas'};
checkFinish   = @(x) any(validatestring(x,validFinishes));
addParameter(p, 'Space', 'sample',checkFinish);
addParameter(p, 'Maxpoints', 8e4, @isnumeric);
parse(p, inputpath, varargin{:});
params = p.Results;
params.Space = validatestring(params.Space, validFinishes);
%==========================================================================
issamp   = false;
switch params.Space
    case 'sample'
        lastbit = "*_locations_sample.mat";
        issamp  = true;
        opts    = loadRegOpts(inputpath);
        if isfield(opts, 'pxsize')
            pxsize = opts.pxsize;
        else
            % slice data: points are in processing pixels in-plane and slice
            % indices along the third column
            pxsize = [opts.processres opts.processres ...
                opts.pxsizes(1)*opts.registres];
        end
        txtuse  = 'Sample space';
    case 'atlas'
        lastbit = "*_locations_atlas.mat";
        txtuse  = 'Atlas space';
end

locpaths = dir(fullfile(inputpath, '**', lastbit));
%==========================================================================
if isempty(locpaths)
    warning("Can't find any valid detections...");
    return
end
%==========================================================================
% cord detections carry the spinal segment of every point, brain ones do not
iscord = ~issamp && ismember('ptsegments', ...
    who('-file', fullfile(locpaths(1).folder, locpaths(1).name)));
if iscord
    cordgrid = makeSpinalCordGrid();
    res      = cordgrid.atlasres;
end
%==========================================================================
Nchans = numel(locpaths);
if iscord
    % a cord is long and thin, so channels go one above the other
    fh = figure('Position',[50, 50, 1400, 100 + 300*Nchans]);
    p  = panel();
    p.pack('v', Nchans);
else
    fh = figure('Position',[50, 50, 500 + 500*(Nchans-1), 700]);
    p  = panel();
    p.pack('h', Nchans);
end
for ipath = 1:numel(locpaths)
    dpcurr = fullfile(locpaths(ipath).folder, locpaths(ipath).name);
    if issamp
        clocs  = load(dpcurr, "cell_locations");
        clocs  = clocs.cell_locations;
    else
        clocs  = load(dpcurr, "atlasptcoords");
        clocs  = clocs.atlasptcoords;
    end
    Ntotal = size(clocs, 1);
    if iscord
        % points beyond the registered stretch of cord are extrapolated and
        % can land far outside the atlas, where they would only stretch the axes
        inatlas = all(clocs(:, 1:3) >= 0.5 & ...
            clocs(:, 1:3) <= cordgrid.atlassize([2 1 3]) + 0.5, 2);
        clocs   = clocs(inatlas, :);
    end
    clocs = subsampplot(clocs(:, 1:3), params.Maxpoints);

    axcurr = p(ipath).select();
    if issamp
        clocs = clocs.*pxsize;
        scatter3(clocs(:,1), clocs(:,2), clocs(:,3), 2, 'filled');
    elseif iscord
        plotSpinalCordGrid(cordgrid, axcurr); hold on;
        scatter3(clocs(:,3)*res(3), clocs(:,1)*res(2), clocs(:,2)*res(1), 2, 'filled');
    else
        plotBrainGrid([], axcurr); hold on;
        scatter3(clocs(:,2), clocs(:,3), clocs(:,1), 2, 'filled');
    end
    axis equal; axis tight;
    ax = gca; ax.ZDir = 'reverse';
    nameuse = strrep(locpaths(ipath).name,'_',' ');

    if iscord
        title(sprintf('%s, %s (%d of %d points inside the atlas)', txtuse, ...
            nameuse, nnz(inatlas), Ntotal))
    else
        title(sprintf('%s, %s', txtuse, nameuse))
    end
end

%==========================================================================
end
