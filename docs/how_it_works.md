# How LightSuite Works

LightSuite maps central-nervous-system samples into a standard atlas and detects and quantifies labelled cells. Whether you start from a whole-brain lightsheet volume, a spinal cord, or a series of wide-field sections, the backbone is the same:

> **downsample → register to the atlas → (optionally) detect cells → export everything in atlas space.**

This page explains the shared ideas once. The [workflow pages](usage_lightsheet_brain.md) give the click-by-click steps for each modality.

---

## The atlas and coordinate space

* Brains are registered to the **Allen Mouse Brain Common Coordinate Framework** (CCF v3, 2020 parcellation; Wang et al., 2020). Spinal cords use the **Fiederling et al. (2021)** atlas.
* All registration runs on volumes downsampled to **20 µm isotropic** — fast to handle, and still finer than the >500 µm spatial scale of the non-rigid fit. The atlas template and annotations are kept at **10 µm**.
* Results are reported in **atlas voxel coordinates**. The exact column order depends on the module and is stated with each output (for example, slice-module cell coordinates are `[ML, AP, DV]`, whereas the probe `probe_ccf` uses `[AP, DV, ML]`).

---

## Registration: from sample to atlas

B-spline deformable registration is non-convex, so a good initialization is critical. Instead of relying on image similarity alone — which can fail on damaged, distorted, or autofluorescent tissue — LightSuite initializes from **sample-derived 3D point clouds** and then refines in three stages:

1. **Point-cloud initialization.** High-saliency structural points are extracted *locally* (so uneven illumination does not bias them) and matched to a point cloud from the atlas with Bayesian Coherent Point Drift. A **similarity** transform (rotation, uniform scale, translation) is fit first — scale matters because clearing can shrink or swell tissue.
2. **Affine.** A global affine transform is fit from matched points (automatic correspondences plus any you add).
3. **B-spline (Elastix).** A smooth, local deformation is computed with a multi-metric objective: image mutual information **plus** a corresponding-points term that pulls your landmarks onto their matches. Both the forward (atlas→sample) and inverse (sample→atlas) transforms are saved; the inverse is what maps detected cells into atlas space.

Each modality adapts this backbone: the spinal cord is **straightened** first, and the slice module registers each 2D section before stacking them into a volume.

---

## Refining registration with control points (active learning)

A shared GUI lets you add corresponding landmarks between your sample and the atlas (the brain and spinal-cord modules use the same 3D interface; the slice module has a 2D version). You click matching points on the sample and atlas side by side, and the atlas region boundaries are overlaid on your sample.

The overlay *is* the feedback loop. Whenever the sample and atlas landmark counts match, a 3D affine is re-fit live from all pairs so far and the overlay is redrawn — so every new point is placed against a progressively better alignment. Candidate views are sampled evenly along all three axes in randomized order (and expose one quadrant at a time) to encourage uniform coverage rather than clustering.

**Tips that apply everywhere:**

* Place **at least 16 pairs** spread across the full extent of the sample; the fit only activates at 16, and 100+ helps on heavily deformed tissue.
* Watch the **MSE** in the banner — lower means a tighter fit.
* Good landmarks are ventricle boundaries, major fiber tracts, and distinct nuclei outlines.
* In damaged or missing regions, place a point where the structure *should* be. This stops the deformation from bending the atlas into the void.

---

## Cell detection and artifact classification

Detection borrows from spike-sorting: extract candidates everywhere, then triage them.

1. **Bandpass filter.** A Butterworth filter centered on the expected cell size highlights cell-sized objects and suppresses noise and slow background.
2. **Signal-to-background ratio (SBR).** Subtracting the bandpass result from the original gives a local background; dividing by it converts raw intensity into an SBR. SBR is robust to the threefold intensity differences seen across brain areas and to one-sided lightsheet illumination.
3. **Local maxima.** Maxima above the primary SBR threshold become candidate centers, grown out to a secondary threshold.
4. **Morphological filtering.** Candidates that are too elongated (axis ratio > 2.5), too small, or too large/bright (tissue folds, bubbles) are discarded.
5. **CNN artifact classifier.** A lightweight network classifies each candidate from three maximum-intensity projections (cheaper than a 3D network) and removes false positives such as bubble edges, tissue–solution interfaces, and bright neurites. It reaches >98% accuracy and can be retrained on your own labels.

The whole-brain pipeline runs this in 3D; the slice module uses a 2D adaptation. `celldiam` sets the filter scale (the single most important parameter) and the threshold pair sets sensitivity.

---

## What LightSuite produces

* **Registration transforms** (`transform_params.mat`) — forward and inverse, reused by every downstream step including implant tracing.
* **Registered volumes** for every channel, warped into atlas space.
* **Cell positions** in atlas space, plus **regional cell counts and fluorescence intensities** per area and hemisphere (exportable to CSV).
* **Implant localization** for optical fibers / GRIN lenses and Neuropixels probes — see [Probe & implant tracing](usage_tracing.md).

---

## References

* Allen CCF v3: Wang et al., *Cell*, 2020.
* Spinal cord atlas: Fiederling et al., *Cell Rep. Methods*, 2021.
* Coherent Point Drift (BCPD): Hirose, [github.com/ohirose/bcpd](https://github.com/ohirose/bcpd).
* Deformable registration: [Elastix](https://elastix.lumc.nl/).
* LightSuite: Karamanlis et al., *high-throughput registration and cell counting in the central nervous system* (in preparation).
