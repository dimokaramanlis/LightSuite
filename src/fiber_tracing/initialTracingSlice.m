function islice = initialTracingSlice(vol, all_points, slice_dim)
%INITIALTRACINGSLICE First slice to show in the tracing annotation GUIs.
%   islice = initialTracingSlice(vol, all_points, slice_dim)
%
%   Returns the median position (along slice_dim) of all annotated points
%   across all groups. If no points exist, returns the median of the slices
%   that contain any non-zero data (relevant for slice datasets, where many
%   AP planes are empty). Falls back to the volume centre.

if nargin < 3; slice_dim = 1; end
nsl = size(vol, slice_dim);

allpts = vertcat(all_points{:});
if ~isempty(allpts)
    islice = round(median(allpts(:, slice_dim)));
else
    other_dims = setdiff(1:ndims(vol), slice_dim);
    hasdata    = squeeze(any(vol, other_dims));
    inz        = find(hasdata);
    if isempty(inz)
        islice = round(nsl / 2);
    else
        islice = round(median(inz));
    end
end

islice = min(max(islice, 1), nsl);
end
