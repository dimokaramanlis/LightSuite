function all_results = fitGRINFibers(savepath, all_points, diameter_um, trstruct, registres, volsize, av)
%FITGRINFIBERS Fit and save GRIN/fiber implants from points clicked on the registration volume.
%
%   ALL_RESULTS = FITGRINFIBERS(SAVEPATH, ALL_POINTS, DIAMETER_UM, TRSTRUCT,
%   REGISTRES, VOLSIZE) maps each group of clicked points to the atlas
%   (registrationPointsToAtlas), fits a cylinder axis per fiber
%   (fitFiberInAtlas), samples the annotation across the axis at
%   grinTracingDepths, and writes grin_fiber<N>_atlas.mat to SAVEPATH. This is
%   what annotateGRINLens does when you press F; it is a function so a fit can
%   be redone from grin_fiber_points.mat without the GUI, e.g. after the
%   registration transform changed.
%
%   ALL_POINTS is the 1 x G cell of N x 3 display-volume points saved by the
%   GUI; groups with fewer than 4 points are skipped. AV (optional) is the
%   annotation volume, loaded if not given. ALL_RESULTS has one struct per
%   fitted fiber, with .group set to its group number.
%
%   See also ANNOTATEGRINLENS, FITFIBERINATLAS, PLOTGRINTRACINGRESULTS.

if nargin < 7 || isempty(av)
    allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
    if isempty(allen_atlas_path)
        error('annotation_10.nii.gz not found on MATLAB path.');
    end
    fprintf('Loading annotation volume...\n');
    av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
end

grp_cnts   = cellfun(@(x) size(x, 1), all_points);
valid_grps = find(grp_cnts >= 4);
all_results = cell(1, numel(valid_grps));
if isempty(valid_grps)
    return;
end
fprintf('--- GRIN Lens Fitting (%d fiber(s)) ---\n', numel(valid_grps));

% Sample the cross-sections above the tip as well as below it (negative
% offsets run back up the fiber track toward the entry). The whole range
% is stored; plotGRINAtlasImages shows only the non-negative depths.
depths_um  = grinTracingDepths();
depths_vox = depths_um / 10;
Ndepths    = numel(depths_um);
radius_vox = (diameter_um / 2) / 10;

% all fibers' points to atlas space in one call (same transform as cells)
counts     = grp_cnts(valid_grps);
atlas_all  = registrationPointsToAtlas(cat(1, all_points{valid_grps}), trstruct, ...
    registres, savepath, volsize);
grp_ends   = cumsum(counts(:)');
grp_starts = [1, grp_ends(1:end-1) + 1];

for ki = 1:numel(valid_grps)
    g = valid_grps(ki);
    fprintf('\n[Fiber %d]  %d points\n', g, counts(ki));
    atlas_pts = atlas_all(grp_starts(ki):grp_ends(ki), :);
    fprintf('  Atlas: AP [%.0f–%.0f]  DV [%.0f–%.0f]  ML [%.0f–%.0f]\n', ...
        min(atlas_pts(:,1)), max(atlas_pts(:,1)), min(atlas_pts(:,2)), max(atlas_pts(:,2)), ...
        min(atlas_pts(:,3)), max(atlas_pts(:,3)));

    % av is passed in so the fit can read the brain surface where the
    % implant enters and keep the axis within 45° of going straight in.
    [center_vox, normal_vox, fit_info] = fitFiberInAtlas(atlas_pts(:, [2 1 3]), radius_vox, av);
    fprintf('  Centre: [%.1f  %.1f  %.1f]  Normal: [%.3f  %.3f  %.3f]\n', center_vox, normal_vox);

    slices_av = cell(1, Ndepths);
    rvec_arr  = cell(1, Ndepths);
    for ii = 1:Ndepths
        [slices_av{ii}, rvec_arr{ii}] = extractAnnotationCircularSlice( ...
            av, center_vox, normal_vox, radius_vox, depths_vox(ii));
    end

    res = struct();
    res.center_vox = center_vox;
    res.normal_vox = normal_vox;
    res.radius_vox = radius_vox;
    res.atlas_pts  = atlas_pts;
    res.depths_um  = depths_um;
    res.slices_av  = slices_av;
    res.rvec_arr   = rvec_arr;
    res.fit_info   = fit_info;

    outfile = fullfile(savepath, sprintf('grin_fiber%d_atlas.mat', g));
    save(outfile, '-struct', 'res');
    fprintf('  Results saved to %s\n', outfile);
    res.group = g;
    all_results{ki} = res;
end
end
