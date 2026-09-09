# Probe & Implant Tracing

Once a brain volume has been registered to the Allen CCF, LightSuite can localize implanted hardware directly in atlas space. Two complementary tools share the same registered volume and the same coordinate transform:

* **Neuropixels (linear) probes** — traced as a straight line of best fit, yielding the ordered list of brain regions along the shank.
* **Cylindrical implants (optical fibers / GRIN lenses)** — traced as a cylinder, yielding the atlas regions beneath the implant tip and the fluorescence intensity below it.

This is the implant-localization output of the LightSuite paper (Fig. 1G and Ext. Data Fig. 1): the implant axis and terminal face in atlas space, together with the estimated field of view at the tip, which resolves directly whether opsin expression falls within the effective photon cone. Both tools operate on the 20 µm registration volume and reuse the registration's [similarity → affine → B-spline transform chain](how_it_works.md#registration-from-sample-to-atlas), so no additional registration is required.

---

## Before You Start

You will need:

* **A fully registered brain volume.** Run the [whole-brain workflow](usage_lightsheet_brain.md) first. The `savepath` folder must contain `regopts.mat`, `transform_params.mat`, and the `chan_X_sample_register_20um.tif` registration volume.
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

**Depth range.** Parcellation labels and fluorescence are sampled on square segments perpendicular to the axis, twice the nominal implant radius across, at `depths_um = -400:100:400` µm (`grinTracingDepths`). 0 is the lens bottom face and positive goes deeper; the **negative** offsets run back up the track toward the entry, a check that the fitted axis follows the visible fiber. Figures draw only the non-negative depths; the negative ones are stored.

Regenerate the figures from disk with `plotGRINTracingResults(savepath)`. The companion intensity figure is produced automatically when a `volume_registered/` folder is present.

That call also **writes the sampled cross-sections back into each `grin_fiber<N>_atlas.mat`**, so they live alongside the fit they came from:

| Field | Contents |
| :--- | :--- |
| `depths_um`, `slices_av`, `rvec_arr` | Atlas cross-sections at the full depth range, re-sampled if the file was saved at a narrower one |
| `slices_int` | 1×Ndepths cell, each `[Ngrid × Ngrid × Nchan]` single — the fluorescence disc at that depth |
| `median_intensity` | `[Ndepths × Nchan]` median inside the implant disc |
| `chan_names` | 1×Nchan channel names, labelling the third dimension of `slices_int` |

The file is rewritten in place, atomically: fit fields are never touched and other variables survive. A file is rewritten only when its depths were re-sampled or intensity was added, so files from an older `annotateGRINLens` (5 depths, no intensity) are upgraded the first time `plotGRINTracingResults` runs on the folder.
### How the axis is fitted

**PCA does not work here.** For a cloud filling a cylinder of length `L` and radius `R`, the variance along the axis is `L²/12` against `R²/4` radially, so the leading eigenvector is the axis only when `L > √3·R` — an implant entering deeper than ~0.87 diameters. Below that PCA returns a *radial* direction, nearly 90° wrong, which is exactly the case of a lens that barely enters the brain. Robust regression does not fix it either: *any* fit minimising summed squared distance to a line **is** PCA (that sum is `trace(S) − nᵀSn`), and a mean-square criterion always prefers the fat direction of a short cloud.

`fitFiberInAtlas` instead takes the trajectory to be **the line minimising the radius of the cylinder enclosing 95% of the clicks**. That minimum sits at the true axis for any aspect ratio — tilting a cylinder inside its covering cylinder can only widen it — and the excluded 5% absorbs misclicks and spill. The width is evaluated as the **mean of the ordered point-to-line distances between the 80th and 95th percentiles**, not a single quantile: a plain quantile lets a tilted axis shed the very end points that reveal the tilt, roughly halving the direction signal, while the averaged band keeps the outer envelope and stays continuous for finite-difference gradients.

The problem is solved by **sequential quadratic programming** (`fmincon`) over four numbers — two tilt angles away from the entry normal and two lateral coordinates of the line — so the ±45° limit is a bound enforced at every iterate rather than clipped afterwards. The objective is non-convex in tilt, so the solve starts from the best points of a coarse scan; `fit_info.at_tilt_bound` records whether the solution ended on the boundary. The implant bottom face is then placed at the **99th percentile** of the axial coordinate of the enclosed points.

### The entry constraint

Implants enter roughly perpendicular to the brain surface, and the fit is told so. The surface is the labelled voxels of the annotation volume that have an unlabelled neighbour; a plane is fit to the patch within **600 µm** of the entry point, and its smallest-variance direction is the normal. Eigen-decomposition is right here and wrong for the axis, because a surface patch is a genuinely flat sheet where a click cloud is not a genuine line. In simulation the normal comes back within **1.5–4°** at mouse-cortex curvature.

*Nearest* surface is not automatically the entry. The interhemispheric fissure and the ventricles are unlabelled sheets *inside* the brain, and for a midline implant their walls are often the closest "surface" — on the real CCF one midline site came back 92° wrong for exactly this reason. A candidate patch must therefore also be **exposed**: walking outward from it must stay in open space. A fissure wall runs back into tissue within a few voxels and is rejected, and the next-nearest patch is tried (up to 8). With that check the same site returns to 5.8°.

The axis may then lean up to `'MaxTilt'` (**45°**) from that normal about each transverse axis. For a long implant the constraint never binds. For a short one ±45° alone is not enough — the fit still wanders tens of degrees on noise — so `'EntryPrior'` adds a **rotation-invariant** quadratic pull back towards the normal, weighted by how much the annotation can say about the lean at all. That weight comes from the ratio of the two largest scatter eigenvalues after trimming strays: 1 for a ring or blob (no direction information), growing with the length of the annotated shaft. It is reported as `fit_info.depth_in_diameters`, the weight as `fit_info.prior_weight` (0 = clicks alone). A blob returns the entry normal; a shaft annotated over a couple of diameters overrules the prior.

> Measuring that depth *along the entry normal* would be self-defeating: a wrong normal foreshortens a deep annotation into an apparently flat one, and the fit then over-weights that same wrong normal on the strength of it.

Passing `av` is what makes this work. Without it (`fitFiberInAtlas(pts, radius)`) the constraint falls back to plain DV, which is only right for a vertical implant on flat cortex.

### Accuracy

Simulated implants (60 clicks through the cylinder, 15 µm jitter, 5% gross misclicks, spherical brain at mouse-cortex curvature). Axis error in degrees, mean over 8 runs:

| annotated depth | lean 0° | lean 10° | lean 25° | PCA fit |
| ---: | ---: | ---: | ---: | ---: |
| 0.4 diameters | 5.1° | 15.0° | 33.0° | 77–81° |
| 0.8 diameters | 5.4° | 9.9° | 22.8° | 66–76° |
| 1.5 diameters | 5.6° | 5.5° | 8.7° | 11–17° |
| 3.0 diameters | 2.1° | 2.2° | 2.6° | 4–5° |

Tip position lands within 35–210 µm; the fit takes ~0.2 s per fiber. On the **real CCF** at five cortical sites, against an independently computed surface normal: 0.3–8.5° at 1.5 diameters, 0.7–3.1° at 3 diameters.

**Depth of annotation buys accuracy, not click count** — though clicks help: at 1.5 diameters, 12 clicks give 18°, 25 give 7.6°, 60 give 4.7°, 150 give 3.3°.

### Tips & Troubleshooting

* **Spread the clicks in depth.** A ring of points at a single depth carries almost no information about the lean; the fit returns the surface normal and warns you. Every extra slice annotated is worth more than extra points on a slice you already did. Click the part of the implant *above* the brain surface too if visible — it extends the annotated depth for free.
* Use at least ~25 points per fiber; below ~12 the axis is poorly determined. Points may sit anywhere on, in or just around the implant, but **do not** click other structures in the same group: a stray beyond the outer 5% is ignored, a systematic second cluster is not.
* **`fit_info.at_tilt_bound`** means the clicks want an implant more oblique than the surface allows. Take it seriously rather than raising `MaxTilt`: usually two implants' clicks in one group, a blob with no direction, or a genuinely oblique insertion near the midline. Check `depth_in_diameters` first — below ~1 the direction is not in the data.
* **`radius_fit`** much larger than your `Diameter` means the clicks spread wider than the implant. The drawn cylinder and the sampled cross-sections always use the **nominal** diameter, so it will look too thin for the cloud.
* `fit_info` also carries `tilt_deg`, `prior_weight`, the inlier mask and the residuals. Warnings are printed at the ±45° limit and when the annotation covers less than one diameter of depth.
* Intensity cross-sections require registered channel volumes in `volume_registered/` (from `generateRegisteredBrainVolumes`).
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
* LightSuite: Karamanlis et al., *Scalable atlas registration and cell detection across the central nervous system* (in preparation).
