# Probe & Implant Tracing

Once a lightsheet brain has been registered to the Allen CCF, LightSuite can localize implanted hardware directly in atlas space. Two complementary tools share the same registered volume and the same coordinate transform:

* **Neuropixels (linear) probes** — traced as a straight line of best fit, yielding the ordered list of brain regions along the shank.
* **Cylindrical implants (optical fibers / GRIN lenses)** — traced as a cylinder, yielding the atlas regions beneath the implant tip and the fluorescence intensity below it.

This corresponds to the implant-localization output in the LightSuite paper (Fig. 1F: *"localization of fiber-like implants in atlas space along with fields of view"*). Both tools operate on the 20 µm registration volume and reuse the registration's [similarity → affine → B-spline transform chain](how_it_works.md#registration-from-sample-to-atlas), so no additional registration is required.

---

## Before You Start

You will need:

* **A fully registered lightsheet brain.** Run the [Lightsheet brain workflow](usage_lightsheet_brain.md) first. The `savepath` folder must contain `regopts.mat`, `transform_params.mat`, and the `chan_X_sample_register_20um.tif` registration volume.
* **Allen CCF atlas resources on the MATLAB path:** `annotation_10.nii.gz` and `parcellation_to_parcellation_term_membership.csv`. The 3-D trajectory plots also require `brainGridData.npy`.

| Tool | Function | Shape fit | Main output |
| :--- | :--- | :--- | :--- |
| Neuropixels probes | `annotateNeuropixelsProbes` | Line (SVD) | `probe_ccf.mat` |
| GRIN lens / optical fiber | `annotateGRINLens` | Cylinder | `grin_fiber<N>_atlas.mat` |

### Choosing the tracing channel

By default both tools display the **registration channel**. Implants are often labelled in a *different* channel (a probe dip-coated in dye, fiber autofluorescence, etc.), so you can trace on any channel:

```matlab
annotateNeuropixelsProbes(savepath, 'Channel', 2);   % trace on channel 2
annotateGRINLens(savepath, 'Channel', 2, 'Diameter', 500);
```

Every channel is downsampled onto the **same 20 µm registration grid**, so the atlas mapping is identical no matter which channel you trace on — only the displayed image changes. Pass `'VolumePath', '<file>'` instead of `'Channel'` to point at a specific volume (absolute path, or a filename inside `savepath`). The active volume is shown in the GUI title bar.

---

## Neuropixels Probe Tracing

**Demo script:** `demos/ls_trace_neuropixels.m`
**Function:** `annotateNeuropixelsProbes(savepath)`

The GUI shows the registration volume in coronal view. You annotate up to nine probes, each as a numbered group, by clicking points along the track wherever it is visible. The track is modelled as a straight line passing through those points, following the approach of [AP_histology](https://github.com/petersaj/AP_histology).

### Workflow

1. **Open the GUI** on a registered folder: `annotateNeuropixelsProbes(savepath)` (add `'Channel', N` to trace on a non-registration channel — see [Choosing the tracing channel](#choosing-the-tracing-channel)).
2. **Select a probe** with number keys `1`–`9` (the active probe is colour-coded).
3. **Click points** along the probe track. Click it in *every* slice where it is visible — points may span many slices.
4. **Press `S`** to save your points (you can reopen and continue later).
5. **Press `F`** to fit. LightSuite transforms all points to Allen CCF space, fits a line of best fit per probe, reads the regions traversed, writes `probe_ccf.mat`, and opens the trajectory figure.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** or **scroll** | Navigate coronal slices |
| **1 – 9** | Select the active probe |
| **Click** | Add a track point to the active probe |
| **Backspace** | Delete the last point of the active probe |
| **c** | Clear the active probe's points on the current slice |
| **Enter** | Jump to a specific slice number |
| **s** | Save points to `neuropixels_probe_points.mat` |
| **f** | Fit all probes, save `probe_ccf.mat`, and plot |

### Outputs

| File | Contents |
| :--- | :--- |
| `neuropixels_probe_points.mat` | Raw clicked points per probe (reloaded automatically on reopen) |
| `probe_ccf.mat` | The fitted `probe_ccf` struct array (see schema below) |

The trajectory figure shows a 3-D brain-grid view with every probe's points and fitted insertion→tip line, plus one region-vs-depth column per probe. You can regenerate it from disk without re-annotating:

```matlab
[probe_ccf, fig] = plotNeuropixelsTracingResults(savepath);
```

### The `probe_ccf` output structure

`probe_ccf` is a struct array with one element per fitted probe. The `points`, `trajectory_coords`, and `trajectory_areas` fields mirror the [AP_histology](https://github.com/petersaj/AP_histology) `probe_ccf` format; all coordinates are Allen CCF **10 µm voxels** in **`[AP, DV, ML]`** order.

| Field | Type | Description |
| :--- | :--- | :--- |
| `probe_number` | scalar | The probe group number (1–9) used in the GUI |
| `points` | N×3 | Clicked points transformed to CCF voxels `[AP, DV, ML]` |
| `trajectory_coords` | 2×3 | `[insertion; tip]` — where the fitted line enters and leaves labelled brain tissue |
| `trajectory_areas` | table | Ordered regions along the trajectory (one row per region) |
| `fit_centroid` | 1×3 | Mean of `points` (line-fit metadata) |
| `fit_direction` | 1×3 | Unit probe direction, oriented ventrally (DV increasing) |
| `fit_endpoints` | 2×3 | Endpoints of the evaluated fit line |

The `trajectory_areas` table has these columns:

| Column | Description |
| :--- | :--- |
| `parcellation_index` | Allen parcellation index of the region |
| `acronym` | Region acronym (e.g. `CA1`, `VISp5`) |
| `name` | Full region name |
| `color` | Region colour, RGB triplet (0–255) |
| `trajectory_depth` | N×2 `[enter, exit]` depth in µm from the brain surface |
| `n_voxels` | Number of sampled voxels in the region span |

> **Compatibility note:** LightSuite labels regions using the Allen 2020 parcellation (`annotation_10.nii.gz`), so `trajectory_areas` carries acronyms/names/RGB colours from that parcellation rather than the 2017 structure-tree rows used by AP_histology. The geometry (`points`, `trajectory_coords`) and the `trajectory_depth` field follow AP_histology exactly.

### Tips & Troubleshooting

* **Click across multiple slices.** A minimum of 2 points defines a line, but spreading several points along the full dorsoventral extent of the track gives a far more stable fit.
* **Aim for the track centre.** The line is fit through your points, so accurate placement matters more than the number of points.
* `Each probe needs at least 2 track points` — add more points to that probe before pressing `F`.
* `probe line never enters a labelled brain region` — your points are likely outside the brain or the registration is off. Re-check the registration overlays before tracing.
* **Re-fitting overwrites `probe_ccf.mat`,** but your clicks live in `neuropixels_probe_points.mat` and are reloaded next time you open the GUI on that folder.
* Probe colours are consistent between the annotation GUI and the trajectory figures.

---

## GRIN Lens / Optical Fiber Tracing

**Function:** `annotateGRINLens(savepath, 'Diameter', 500)`

This tool localizes a cylindrical implant and reports the atlas regions beneath it. After fitting, it shows the parcellation regions on circular cross-sections at a range of depths below the lens face (0–400 µm), and — when registered channel volumes are available — the fluorescence intensity beneath the implant (the *"quantification of intensity below the implant"* described in the paper).

### Workflow

1. `annotateGRINLens(savepath, 'Diameter', 500)` — set `Diameter` (µm) to your fiber/lens diameter (default 500). Add `'Channel', N` to trace on a non-registration channel (see [Choosing the tracing channel](#choosing-the-tracing-channel)).
2. Select a fiber with `1`–`9` and click points marking where the implant is — **on its edge, inside it, just around it** — on every slice where you can see it. The points do not have to lie on the outline; what matters is that they cover the implant over as much **depth** as it is visible.
3. Press `S` to save, then `F` to fit. Each fiber is fit as a cylinder in CCF space and its sub-lens atlas regions are displayed.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** or **scroll** | Navigate coronal slices |
| **1 – 9** | Select the active fiber |
| **Click** | Add a point on the active fiber |
| **Backspace** | Delete the last point of the active fiber |
| **c** | Clear the active fiber's points on the current slice |
| **Enter** | Jump to a specific slice number |
| **s** | Save points to `grin_fiber_points.mat` |
| **f** | Fit all fibers and show atlas (and intensity) figures |

### Outputs

| File | Contents |
| :--- | :--- |
| `grin_fiber_points.mat` | Raw clicked points per fiber (`all_points`) |
| `grin_fiber<N>_atlas.mat` | Per-fiber fit: `center_vox`, `normal_vox`, `radius_vox`, `atlas_pts`, `depths_um`, `slices_av`, `rvec_arr`, `fit_info` |

**Depth range.** Cross-sections are sampled at `depths_um = -400:100:400` µm along the fitted axis (`grinTracingDepths`). 0 is the lens bottom face, positive goes deeper into the brain, and the **negative** offsets go back up the fiber track toward the entry — sampling above the tip is a check that the fitted axis really follows the visible fiber. The **figures show only the non-negative depths** (0–400 µm), so they look the same as before; the negative depths are stored, not drawn.

Regenerate the figures from disk with `plotGRINTracingResults(savepath)`. The companion intensity figure is produced automatically when a `volume_registered/` folder is present.

That call also **writes the sampled cross-sections back into each `grin_fiber<N>_atlas.mat`**, so they live alongside the fit they came from:

| Field | Contents |
| :--- | :--- |
| `depths_um`, `slices_av`, `rvec_arr` | Atlas cross-sections at the full depth range, re-sampled if the file was saved at a narrower one |
| `slices_int` | 1×Ndepths cell, each `[Ngrid × Ngrid × Nchan]` single — the fluorescence disc at that depth |
| `median_intensity` | `[Ndepths × Nchan]` median inside the implant disc |
| `chan_names` | 1×Nchan channel names, labelling the third dimension of `slices_int` |

The file is rewritten in place (atomically), the fit fields are never touched, and any other variables survive the rewrite. A file is only rewritten when its depths were re-sampled or intensity was added, and re-running refreshes these fields. Files saved by an older `annotateGRINLens` (5 depths, no intensity) are upgraded to the full range the first time `plotGRINTracingResults` runs on the folder.

### How the axis is fitted

**A PCA fit is not usable here.** For a cloud filling a cylinder of length `L` and radius `R`, the variance along the axis is `L²/12` against `R²/4` in each radial direction, so the leading eigenvector is the axis only when `L > √3·R`, i.e. when the implant enters deeper than ~0.87 diameters. Below that PCA returns a *radial* direction and the axis comes out nearly 90° wrong — exactly the case of a GRIN lens that barely enters the brain.

**Nor would a robust regression fix it,** because the failure is not caused by outliers. *Any* fit minimising the summed squared distance to a line **is** PCA — that sum equals `trace(S) − nᵀSn`, so minimising it maximises `nᵀSn` — and a mean-square criterion always prefers the fat direction of a short cloud. What is needed is a criterion measuring how **wide** the cloud is around the line, not how far its points are on average.

`fitFiberInAtlas` therefore takes the axis to be **the line minimising the radius of the cylinder that covers 95% of the clicks**. That minimum sits at the true axis for any aspect ratio — tilting a cylinder inside its covering cylinder can only widen it — and it is robust by construction, since the outermost 5% (misclicks, spill around the implant) are excluded.

This is posed as a **constrained optimisation in four numbers** and solved with `fmincon`: the two tilt angles of the direction away from the entry normal, and the two coordinates of the line's lateral position. The ±45° limit is then a plain bound that the solver enforces at every iterate rather than something clipped afterwards. The objective is not convex in the tilt — a cloud can look narrow around more than one direction — so the solve starts from the best few points of a coarse scan of the tilt box. `fit_info.at_tilt_bound` records whether the solution ended up on the boundary.

One detail matters more than it looks: the width is not a single quantile but the **mean of the outermost distances up to that quantile**. Everything above a plain quantile is discarded, so a tilted axis can simply *shed* the end points that give the tilt away, and the direction signal roughly halves. Averaging the band just below it keeps the outer envelope — where the direction information lives — while the 5% above absorbs the strays. It also makes the objective an L-statistic, so it is continuous (the two clicks that swap rank are equidistant at the swap) and finite-difference gradients behave.

### The entry constraint

An implant is inserted roughly perpendicular to the brain surface, and the fit is told so. The annotation volume is `0` outside the brain and `>0` inside, so the surface is simply the labelled voxels having an unlabelled neighbour. The patch of it nearest the clicks is locally flat, and the normal of a plane *is* its smallest-variance direction — so an eigen-decomposition is the right tool for the surface even though it is the wrong tool for the axis, a surface patch being a genuinely flat sheet where the click cloud is not a genuine line. In simulation the entry normal comes back **within 1.5–4° of truth** at mouse-cortex curvature.

*Nearest* surface is not automatically the entry, though. The interhemispheric fissure and the ventricles are unlabelled sheets inside the brain, and their walls are "surface" too — for an implant near the midline they are often the **closest** surface, and fitting a plane to them returns the wall normal, about 90° from the right answer. On the real CCF a midline site came back 92° wrong for exactly this reason. So a candidate patch also has to be **exposed**: walking outward from it must stay in open space, which is where the implant came from. A fissure wall runs back into tissue within a few voxels and is rejected, and the next-nearest patch is tried (up to 8). With that check the same midline site returns to 5.8°.

The axis may then lean up to **`'MaxTilt'` (45°) from that normal in each of the two transverse axes**. For a long implant the constraint never binds and the clicks decide. For a short one, ±45° on its own is *not* enough — inside that box the fit still wanders tens of degrees on noise — so `'EntryPrior'` adds a quadratic pull back towards the entry normal, weighted by how much the annotation can say about the lean at all.

That weight comes from the **shape of the click cloud, measured rotation-invariantly**: the ratio of the two largest scatter eigenvalues is 1 for a ring or a blob of clicks, which carries no direction information, and grows with the length of the annotated shaft. Strays are trimmed first — a covariance is not robust, and three misclicks at arm's length halve the apparent elongation. It is reported as `fit_info.depth_in_diameters`, and the resulting weight as `fit_info.prior_weight` (0 = clicks alone). A blob returns the entry normal; a shaft annotated over a couple of diameters overrules the prior outright.

Measuring that depth *along the entry normal* instead is a trap worth naming, because it is self-defeating: a wrong entry normal foreshortens a deep annotation into an apparently flat one, and the fit then over-weights that same wrong normal on the strength of it.

Passing `av` is what makes all of this work. Without it (`fitFiberInAtlas(pts, radius)`) the constraint falls back to plain DV, which is only right for a vertical implant on a flat piece of cortex.

### Accuracy

Simulated implants, 60 clicks scattered through the cylinder with 15 µm jitter and 5% gross misclicks, on a spherical brain of mouse-cortex curvature. Axis error in degrees, mean over 8 runs:

| annotated depth | lean 0° | lean 10° | lean 25° | old PCA fit |
| ---: | ---: | ---: | ---: | ---: |
| 0.4 diameters | 5.1° | 15.0° | 33.0° | 77–81° |
| 0.8 diameters | 5.4° | 9.9° | 22.8° | 66–76° |
| 1.5 diameters | 5.6° | 5.5° | 8.7° | 11–17° |
| 3.0 diameters | 2.1° | 2.2° | 2.6° | 4–5° |

Tip position lands within 35–210 µm. The fit takes ~0.2 s per fiber. **Depth of annotation is what buys accuracy, not the number of clicks** — though clicks help too: at 1.5 diameters, 12 clicks give 18°, 25 give 7.6°, 60 give 4.7° and 150 give 3.3°.

Repeated on the **real CCF** at five cortical sites (midline, lateral, anterior, posterior), against an independently computed surface normal, the same pattern holds: 0.3–8.5° at 1.5 diameters and 0.7–3.1° at 3 diameters.

### When the fit is fighting you

`fit_info.at_tilt_bound` means the solution sits on the ±45° boundary: the clicks want an implant more oblique than the brain surface allows. That is worth taking seriously rather than raising `MaxTilt`, because in practice it means one of three things — the clicks belong to more than one structure, the annotation is a blob with no real direction in it, or the implant genuinely did not go in perpendicular to the surface (common near the midline, where the surface slopes steeply into the interhemispheric fissure while the implant went in vertically).

Check `depth_in_diameters` first. Below ~1, the direction is not really in the data and the fit is mostly reporting the brain surface. Also compare `radius_fit` against your `Diameter`: much larger means the clicks are spread wider than the implant, which both weakens the direction and makes the drawn cylinder — which is always the **nominal** size — look too thin for the point cloud in the figures.

`fit_info` carries `tilt_deg` (the per-axis lean from the entry normal), `at_tilt_bound`, `depth_in_diameters`, `prior_weight`, `radius_fit` (the covering radius of the click cloud — reported only; the cross-sections are still sampled at the nominal `Diameter`), the inlier mask and the residuals. Two warnings are printed: one when the axis runs into the ±45° limit (the clicks want an implant far more oblique than the surface allows — usually two implants' clicks in one group), and one when the annotation covers less than one diameter of depth.

### Tips

* **Spread the clicks in depth — that is what determines the axis.** A ring of points at a single depth carries almost no information about the lean, and the fit will simply return the brain-surface normal and warn you. Every extra slice you annotate is worth more than extra points on a slice you already did.
* **Click the part of the implant above the brain surface too, if it is visible.** It extends the annotated depth at no cost, which is exactly the quantity that pins the axis down.
* Points may sit anywhere on, in or just around the implant — the fit does not assume they lie on the outline. **Do not** click other structures in the same group: a stray beyond the outer 5% is ignored, but a systematic second cluster will tilt the fit.
* Use at least ~25 points per fiber; below ~12 the axis is poorly determined.
* The intensity cross-sections require registered channel volumes in `volume_registered/` (produced by `generateRegisteredBrainVolumes`).

---

## Output File Summary

```
<savepath>/
├── neuropixels_probe_points.mat   # Raw probe clicks (Neuropixels)
├── probe_ccf.mat                  # Fitted probe_ccf struct array
├── grin_fiber_points.mat          # Raw fiber clicks (GRIN)
└── grin_fiber<N>_atlas.mat        # Per-fiber cylinder fit + sub-lens regions
```

---

## References

* Probe `probe_ccf` format and line-fit convention: **AP_histology**, [github.com/petersaj/AP_histology](https://github.com/petersaj/AP_histology).
* Atlas coordinate framework: Allen Mouse Brain Common Coordinate Framework (Wang et al., *Cell*, 2020).
* LightSuite: Karamanlis et al., *high-throughput registration and cell counting in the central nervous system* (in preparation).
