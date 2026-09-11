function griddata = makeSpinalCordGrid(varargin)
%MAKESPINALCORDGRID Wire mesh of the spinal cord atlas, the cord brainGridData.
%
%   GRIDDATA = MAKESPINALCORDGRID() traces a wire mesh of the spinal cord atlas
%   on the MATLAB path, for PLOTSPINALCORDGRID to draw points against. There is
%   no ready-made mesh for the cord the way brainGridData is for the brain, so it
%   is built from the annotation volume:
%
%     rings      - the outline of the cord, and of its gray matter, in the
%                  transverse plane at the centre of every spinal segment and at
%                  both ends of the cord. The atlas holds one reference section
%                  per segment, so these planes carry all of its cross-sections.
%     long lines - NLONG lines along the cord, joining the same angular position
%                  on consecutive rings, which turn the rings into a tube.
%
%   GRIDDATA is a struct with fields
%     outline    - N x 3 [x y z] in micrometres, rings and long lines one after
%                  the other, separated by NaN rows
%     greymatter - M x 3, the same for the gray matter outline of every ring
%     segments   - table with the name of every spinal segment and, in um along
%                  z, where it starts, ends and where its ring sits
%     atlasres   - native voxel size, [10 10 20] um per [row col plane]
%     atlassize  - native volume size, [rows cols planes]
%
%   The mesh is in micrometres of the native atlas grid, [x y z] = [col row
%   plane] .* ATLASRES([2 1 3]), so atlas-space points as
%   TRANSFORMCORDPOINTSTOATLAS saves them sit on it after the same scaling.
%
%   Reading the atlas takes a few seconds, so the grid is kept for the rest of
%   the MATLAB session.
%
%   Optional name-value arguments
%     nlong  - number of long lines (default 8)
%     reload - rebuild even when a grid is cached (default false)
%
%   See also PLOTSPINALCORDGRID, PLOTBRAINGRID, LOADSPINALCORDATLAS.

p = inputParser;
addParameter(p, 'nlong', 8, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'reload', false, @(x) islogical(x) || isscalar(x));
parse(p, varargin{:});
nlong = round(p.Results.nlong);

persistent cached cachednlong
if ~p.Results.reload && ~isempty(cached) && isequal(cachednlong, nlong)
    griddata = cached;
    return
end
%--------------------------------------------------------------------------
[~, av, parcelinfo, segmentinfo, atlasres, nativesize] = loadSpinalCordAtlas();

% gray matter is whatever the curated hierarchy puts under GM; the central
% canal sits outside it, which is expected and not worth a warning here
wstate = warning('off', 'cordAreaHierarchy:unassignedRegions');
hier   = cordAreaHierarchy(parcelinfo, unique(av(:)));
warning(wstate);
gmids  = hier.parcellation_index(strcmp(hier.division_acronym, 'GM'));
%--------------------------------------------------------------------------
% one ring per segment centre, plus the two ends so the tube is closed
segstart = double(segmentinfo.Start(:));
segend   = double(segmentinfo.End(:));
segring  = round((segstart + segend) / 2);
Nz       = size(av, 3);
ringplanes = unique([1; segring; Nz]);
Nrings     = numel(ringplanes);

if ismember('Segment', segmentinfo.Properties.VariableNames)
    segnames = string(segmentinfo.Segment);
else
    segnames = string((1:numel(segstart))');
end

% [row col] of plane iz -> [x y z] um; plane p spans (p +- 0.5)*res
toum = @(rc, iz) [rc(:, 2)*atlasres(2), rc(:, 1)*atlasres(1), ...
    iz*atlasres(3)*ones(size(rc, 1), 1)];

thetas  = 2*pi*(0:nlong-1)'/nlong;
anchors = nan(nlong, 3, Nrings);
rings   = cell(Nrings, 1);
gmrings = cell(Nrings, 1);
for ir = 1:Nrings
    iz    = ringplanes(ir);
    slice = av(:, :, iz);

    % the cord is one piece; keep only its largest outline
    cordloops = traceOutlines(slice > 0);
    [~, imax] = max(cellfun(@(x) size(x, 1), cordloops));
    outer     = cordloops{imax};
    rings{ir} = [toum(outer, iz); nan(1, 3)];

    % gray matter can come apart into islands caudally, so keep every piece
    gmloops     = traceOutlines(ismember(slice, gmids));
    gmloops     = cellfun(@(x) [toum(x, iz); nan(1, 3)], gmloops, 'UniformOutput', false);
    gmrings{ir} = cat(1, zeros(0, 3), gmloops{:});

    % where each long line crosses this ring: the outline point closest in
    % angle, seen from the middle of the cord
    [rr, cc] = find(imfill(slice > 0, 'holes'));
    angs     = atan2(outer(:, 1) - mean(rr), outer(:, 2) - mean(cc));
    for il = 1:nlong
        [~, ipt] = min(abs(angle(exp(1i*(angs - thetas(il))))));
        anchors(il, :, ir) = toum(outer(ipt, :), iz);
    end
end

longlines = cell(nlong, 1);
for il = 1:nlong
    longlines{il} = [reshape(anchors(il, :, :), 3, [])'; nan(1, 3)];
end
%--------------------------------------------------------------------------
griddata            = struct();
griddata.outline    = cat(1, rings{:}, longlines{:});
griddata.greymatter = cat(1, gmrings{:});
griddata.segments   = table(segnames, (segstart - 0.5)*atlasres(3), ...
    (segend + 0.5)*atlasres(3), segring*atlasres(3), ...
    'VariableNames', {'name', 'zstart_um', 'zend_um', 'zring_um'});
griddata.atlasres   = atlasres;
griddata.atlassize  = nativesize;

cached      = griddata;
cachednlong = nlong;
end

%==========================================================================
function loops = traceOutlines(mask)
%TRACEOUTLINES Smoothed closed outlines, [row col], of the pieces of a mask.
%   Holes are filled first - the central canal would otherwise add a ring of
%   its own - and specks too small to read as a shape are dropped.
minpts = 20;
bounds = bwboundaries(imfill(mask, 'holes'), 'noholes');
bounds = bounds(cellfun(@(x) size(x, 1), bounds) >= minpts);
loops  = cellfun(@smoothLoop, bounds, 'UniformOutput', false);
end

%--------------------------------------------------------------------------
function rc = smoothLoop(rc)
%SMOOTHLOOP Take the pixel staircase out of a closed outline.
rc = rc(1:end-1, :);                  % bwboundaries repeats the first point
n  = size(rc, 1);
w  = min(5, n);
rc = movmean([rc(end-w+1:end, :); rc; rc(1:w, :)], w, 1);   % wrap around
rc = rc(w+1:w+n, :);
rc = rc(1:2:end, :);                  % every other point is plenty at 10 um
rc = [rc; rc(1, :)];                  % closed again
end
