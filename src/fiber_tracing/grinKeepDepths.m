function results = grinKeepDepths(results, keepfun)
% GRINKEEPDEPTHS  Restrict the per-depth fields of GRIN result structs.
%
%   RESULTS = grinKeepDepths(RESULTS, KEEPFUN)
%
%   RESULTS  – 1×Nfibers cell array of fiber result structs (atlas or
%              intensity), each with a .depths_um field.
%   KEEPFUN  – function handle applied to depths_um returning a logical mask,
%              e.g. @(d) d >= 0 to drop the negative (above-the-tip) depths.
%
%   Every field indexed by depth is subset to the kept depths: the cell arrays
%   slices_av, rvec_arr and slices_int, and the rows of the [Ndepths×Nchan]
%   matrices median_intensity and mean_intensity.  All other fields (the fit, atlas_pts, ...) are
%   passed through untouched.  Used by the plotting functions to show only the
%   depths that should appear in the figures while the saved files keep the
%   full range.

    for k = 1:numel(results)
        r = results{k};
        if ~isfield(r, 'depths_um') || isempty(r.depths_um); continue; end

        keep = logical(keepfun(r.depths_um));
        nd   = numel(r.depths_um);
        r.depths_um = r.depths_um(keep);

        for f = {'slices_av', 'rvec_arr', 'slices_int'}
            if isfield(r, f{1}) && iscell(r.(f{1})) && numel(r.(f{1})) == nd
                r.(f{1}) = r.(f{1})(keep);
            end
        end
        for f = {'median_intensity', 'mean_intensity'}
            if isfield(r, f{1}) && size(r.(f{1}), 1) == nd
                r.(f{1}) = r.(f{1})(keep, :);
            end
        end

        results{k} = r;
    end
end
