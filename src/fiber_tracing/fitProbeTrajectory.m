function probe = fitProbeTrajectory(points, atlas)
% FITPROBETRAJECTORY  Fit a straight Neuropixels probe trajectory to a set of
%   atlas-space points and read off the brain regions it passes through.
%
%   PROBE = fitProbeTrajectory(POINTS, ATLAS)
%
%   Implements the AP_histology probe-fitting algorithm (Peters et al.,
%   https://github.com/petersaj/AP_histology): a line of best fit is obtained
%   from the singular value decomposition of the mean-centred points, the
%   Allen CCF annotation volume is sampled at ~1 µm steps along that line, and
%   contiguous runs of the same parcellation region are collapsed into an
%   ordered list of areas with their depth extents.
%
%   Inputs
%   ------
%   POINTS – N×3 probe points in Allen CCF 10 µm atlas voxels, column order
%            [AP DV ML].  At least 2 points are required.
%   ATLAS  – struct from loadTracingAtlas (fields av, areaidx, acronyms,
%            names, colors).
%
%   Output struct PROBE (AP_histology "probe_ccf" entry plus fit metadata)
%   ----------------------------------------------------------------------
%   .points            – the input POINTS ([AP DV ML]).
%   .trajectory_coords – 2×3 [entry; tip] CCF coordinates where the fitted
%                        line enters and leaves labelled brain tissue.
%   .trajectory_areas  – table of regions along the trajectory with columns
%                        parcellation_index, acronym, name, color (RGB 0–255),
%                        trajectory_depth (N×2 [enter exit] µm from surface)
%                        and n_voxels.
%   .fit_centroid      – 1×3 mean of POINTS ([AP DV ML]).
%   .fit_direction     – 1×3 unit probe direction ([AP DV ML]), oriented so DV
%                        increases (pointing ventrally / into the brain).
%   .fit_endpoints     – 2×3 endpoints of the evaluated fit line ([AP DV ML]).

    if size(points, 1) < 2
        error('fitProbeTrajectory: need at least 2 points to fit a probe line.');
    end

    av = atlas.av;
    sz = size(av);

    %------------------------------------------------------------------
    % 1. Line of best fit through the points (primary SVD axis)
    %------------------------------------------------------------------
    r0  = mean(points, 1);
    xyz = points - r0;
    [~, ~, V] = svd(xyz, 0);
    direction = V(:, 1)';

    % Ensure the probe points down (DV is column 2; increasing DV = ventral)
    if direction(2) < 0
        direction = -direction;
    end

    % Evaluate the line well beyond the brain in both directions
    line_eval      = [-1000; 1000];
    probe_fit_line = line_eval .* direction + r0;     % 2×3 endpoints

    %------------------------------------------------------------------
    % 2. Sample the annotation volume along the line (~1 µm = 0.1 vox steps)
    %------------------------------------------------------------------
    seg_len_vox = norm(diff(probe_fit_line, 1, 1));   % length in atlas voxels
    n_coords    = max(2, round(seg_len_vox * 10));

    ap = round(linspace(probe_fit_line(1,1), probe_fit_line(2,1), n_coords));
    dv = round(linspace(probe_fit_line(1,2), probe_fit_line(2,2), n_coords));
    ml = round(linspace(probe_fit_line(1,3), probe_fit_line(2,3), n_coords));

    in_bounds = ap >= 1 & ap <= sz(1) & ...
                dv >= 1 & dv <= sz(2) & ...
                ml >= 1 & ml <= sz(3);
    ap = ap(in_bounds); dv = dv(in_bounds); ml = ml(in_bounds);
    if isempty(ap)
        error('fitProbeTrajectory: probe line does not intersect the atlas volume.');
    end

    coords        = [ap(:), dv(:), ml(:)];
    idx           = sub2ind(sz, ap(:), dv(:), ml(:));
    areas_sampled = double(av(idx));

    %------------------------------------------------------------------
    % 3. Collapse contiguous runs into ordered regions (boundaries)
    %------------------------------------------------------------------
    bins       = [1; find(diff(areas_sampled) ~= 0) + 1; numel(idx)];
    boundaries = [bins(1:end-1), bins(2:end)];        % [run_start, run_end] sample idx
    run_area   = areas_sampled(boundaries(:, 1));

    store = run_area > 0;                             % keep in-brain regions only
    if ~any(store)
        error('fitProbeTrajectory: probe line never enters a labelled brain region.');
    end

    %------------------------------------------------------------------
    % 4. Depths along the trajectory (µm), referenced to the brain entry
    %------------------------------------------------------------------
    step_um   = (seg_len_vox * 10) / (n_coords - 1);  % atlas voxel = 10 µm
    first_in  = find(store, 1, 'first');
    last_in   = find(store, 1, 'last');
    depth_ref = boundaries(first_in, 1);

    run_start = boundaries(store, 1);
    run_end   = boundaries(store, 2);
    % trajectory_depth: N×2 [enter, exit] depth per region (µm from brain
    % surface), matching the AP_histology trajectory_areas.trajectory_depth field
    trajectory_depth = ([run_start, run_end] - depth_ref) * step_um;
    n_voxels         = run_end - run_start + 1;
    ids              = run_area(store);

    %------------------------------------------------------------------
    % 5. Map parcellation indices to acronyms / names / colours
    %------------------------------------------------------------------
    nA       = numel(ids);
    acronym  = cell(nA, 1);
    name     = cell(nA, 1);
    color    = zeros(nA, 3);
    for ii = 1:nA
        m = find(atlas.areaidx == ids(ii), 1);
        if ~isempty(m)
            acronym{ii} = atlas.acronyms{m};
            name{ii}    = atlas.names{m};
            color(ii,:) = atlas.colors(m, :);
        else
            acronym{ii} = num2str(ids(ii));
            name{ii}    = '';
            color(ii,:) = [0 0 0];
        end
    end

    parcellation_index = ids;
    trajectory_areas = table(parcellation_index, acronym, name, color, ...
        trajectory_depth, n_voxels);

    %------------------------------------------------------------------
    % 6. Package output (AP_histology probe_ccf fields + fit metadata)
    %------------------------------------------------------------------
    probe.points            = points;
    probe.trajectory_coords = coords([boundaries(first_in,1); boundaries(last_in,2)], :);
    probe.trajectory_areas  = trajectory_areas;
    probe.fit_centroid      = r0;
    probe.fit_direction     = direction;
    probe.fit_endpoints     = probe_fit_line;
end
