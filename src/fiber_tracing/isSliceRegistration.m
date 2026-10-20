function tf = isSliceRegistration(trstruct)
% ISSLICEREGISTRATION  True if TRSTRUCT (transform_params.mat) was written by
%   the slice module (registerSlicesToAtlas) rather than the light-sheet
%   whole-brain pipeline.
%
%   Slice-registered brains have per-slice 2-D affines and no
%   how_to_perm/ori_pxsize; their volume_registered/ is already in atlas
%   space, so tracing runs on it directly (see loadTracingRegistration).

    tf = isfield(trstruct, 'tformaffine_tform_atlas_to_image') && ...
         isfield(trstruct, 'sliceids_in_sample_space');
end
