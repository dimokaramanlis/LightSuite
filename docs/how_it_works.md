# How LightSuite Works

LightSuite performs **atlas registration and cell detection across the central nervous system**. Whether you start from a whole-brain volume, a spinal cord, a series of wide-field sections, or an *in vivo* fUS recording, the backbone is the same:

> **downsample → register to the atlas → (optionally) detect cells → export everything in atlas space.**

Two ideas carry the pipeline: **damage-resilient point-cloud initialization** and **user-steered non-rigid registration** for alignment, and **local signal-to-background normalization** for detection. This page explains them once; the [workflow pages](usage_lightsheet_brain.md) give the click-by-click steps.

---

## The atlas and coordinate space

* Brains are registered to the **Allen Mouse Brain Common Coordinate Framework** (CCFv3; Wang et al., 2020). Spinal cords use the **Fiederling et al. (2021)** atlas; fUS volumes use a vascular template (Brunner et al., 2021) that has itself been registered into CCFv3.
* Registration runs on volumes downsampled to **20 µm isotropic** — far finer than the several-hundred-µm scale of the non-rigid fit, which matches how much clearing and fixation actually deform tissue. Atlas template and annotations are kept at **10 µm**.
* Results are reported in **atlas voxel coordinates**. Column order depends on the module and is stated with each output (slice-module cell coordinates are `[ML, AP, DV]`; the probe `probe_ccf` uses `[AP, DV, ML]`).
* Quantification defaults to the **substructure level** of the Allen parcellation, per hemisphere.

---

## Registration: from sample to atlas

B-spline deformable registration is non-convex, so initialization decides the outcome. Image similarity alone fails on damaged, distorted, or autofluorescent tissue, so LightSuite starts from **point clouds** and refines in three stages:

1. **Point-cloud initialization.** High-saliency structural points are extracted *locally* (so uneven illumination does not bias them) and matched against an atlas cloud with **Bayesian Coherent Point Drift**. A **similarity** transform (rotation, uniform scale, translation) is fit first — scale matters because clearing shrinks or swells tissue. Tears, compression and bubble voids become statistical outliers, and missing structures (a detached olfactory bulb, say) simply remove points instead of dragging the rest of the anatomy with them.
2. **Affine.** A global affine transform is fit by least squares from the matched point set — automatic correspondences plus any landmarks you added.
3. **B-spline (Elastix 5.1.0).** A free-form deformation on a coarse control-point grid (**0.64 mm** for mouse brains) is optimized against a **multi-metric objective**: Advanced Mattes mutual information (48 histogram bins, weight 1.0) **plus** a corresponding-points Euclidean-distance term (weight 0.2) that pulls your landmarks onto their matches. The coarse grid limits free parameters and prevents overfitting without a displacement penalty.

Forward (atlas→sample) and inverse (sample→atlas) transforms are both saved; the inverse maps detected cells into atlas space.

Each modality adapts this backbone: the spinal cord is **straightened and untwisted** first, the slice module registers each 2D section before stacking, and fUS sessions are first aligned to a within-mouse seed.

---

## Refining registration with control points (active learning)

A shared GUI lets you add corresponding landmarks between sample and atlas (brain and spinal cord use the same 3D interface; the slice module has a 2D version). You click matching points side by side, with atlas region boundaries overlaid on your sample.

The overlay *is* the feedback loop. Whenever the sample and atlas landmark counts match, a 3D affine is re-fit live from every pair so far and the boundaries are redrawn — so each new point is placed against a better alignment, and residual mismatches become easier to see. Candidate views are sampled evenly along all three axes in randomized order, one quadrant at a time, to spread coverage rather than let it cluster.

This is what image similarity cannot supply on its own: as landmarks accumulate, the landmark residual keeps falling while mutual information stays flat — so mutual information alone cannot tell an anatomically faithful fit from a plausible but wrong one, and the landmark residual is the more sensitive quality measure.

**Tips that apply everywhere:**

* The fit activates at **16 pairs**; the landmark error keeps improving up to about **100**, which takes 10–20 minutes per sample.
* Watch the **MSE** in the banner — lower means a tighter fit.
* Good landmarks are ventricle boundaries, major fiber tracts, and distinct nuclei outlines — where most annotation naturally concentrates.
* In damaged or missing regions, place a point where the structure *should* be. This stops the deformation from bending the atlas into the void.

---

## Cell detection and artifact classification

Detection borrows from spike-sorting: extract candidates everywhere, then triage them. It runs in GPU-accelerated 3D batches, so a whole brain fits on a standard workstation without global illumination flat-fielding.

1. **Band-pass filter.** A 3D Butterworth filter centered on the expected cell size isolates cell-sized features; the surrounding spatial frequencies estimate the local background.
2. **Signal-to-background ratio (SBR).** Dividing the band-pass volume by that background converts raw intensity into an SBR, factoring out the up-to-threefold baseline differences across brain areas and any illumination asymmetry, such as that of a one-sided light sheet.
3. **Local maxima.** Maxima above the primary SBR threshold become candidate centers, grown out to a secondary threshold.
4. **Morphological filtering.** Candidates too elongated (principal-axis ratio > 2.5), too small, or too large and bright (tissue folds, bubbles) are discarded.
5. **CNN artifact classifier.** A lightweight network classifies each candidate from three orthogonal maximum-intensity projections (XY, XZ, YZ), avoiding the memory cost of 3D convolutions. It reaches **98.7%** validation accuracy and discards **10–20%** of candidates (median 16%) — bubble edges, tissue–solution interfaces, neurite varicosities — which are systematically dimmer, smaller and more elongated than real somata. See [CNN cell classification](usage_lightsheet_brain.md#7-cnn-cell-classification) to label, train and apply your own.

The pipeline sees a stitched 3D stack rather than a microscope, so it applies to any volumetric modality — light-sheet, fMOST, serial two-photon tomography — with `celldiam` and `pxsize` as the only modality-specific settings. The whole-brain pipeline runs this in 3D; the slice module uses a 2D adaptation. `celldiam` sets the filter scale (the single most important parameter) and the threshold pair sets sensitivity.

Counting cells rather than integrating intensity also removes illumination artifacts: with one-sided light sheets, regional intensities are markedly asymmetric across hemispheres while cell counts are not.

---

## What LightSuite produces

* **Registration transforms** (`transform_params.mat`) — forward and inverse, reused by every downstream step including implant tracing.
* **Registered volumes** for every channel, warped into atlas space.
* **Cell positions** in atlas space, plus **regional cell counts and fluorescence intensities** per area and hemisphere (exportable to CSV).
* **Implant localization** for optical fibers / GRIN lenses and Neuropixels probes — see [Probe & implant tracing](usage_tracing.md).

---

## Performance

Measured on a workstation with an Intel Core i9-10900X, 64 GB RAM and an NVIDIA RTX A4000 (16 GB VRAM):

| Stage | Time |
| :--- | :--- |
| Whole-brain cell detection | 3–4 h |
| Automated point-cloud registration | 5–6 min |
| Interactive landmark curation | 10–20 min (~100 landmarks) |

A typical mouse brain yields ~4.5 × 10⁵ detected cells.
