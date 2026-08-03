function depths_um = grinTracingDepths()
% GRINTRACINGDEPTHS  Standard depth offsets (µm) for GRIN / fiber cross-sections.
%
%   DEPTHS_UM = grinTracingDepths()
%
%   Offsets along the fitted fiber axis at which the atlas parcellation and the
%   channel fluorescence are sampled and stored.  0 is the lens bottom face;
%   positive goes deeper into the brain, negative goes back up the fiber track
%   toward the entry.
%
%   The negative offsets are stored so the fit can be checked against the
%   visible fiber above the tip — the cross-sections there should still land on
%   the track if the axis is right.  The plotting functions
%   (plotGRINAtlasImages, plotGRINIntensityImages) show only the non-negative
%   depths, so the figures are unchanged; the negative ones live in the saved
%   grin_fiber<N>_atlas.mat only.

    depths_um = -400 : 100 : 400;
end
