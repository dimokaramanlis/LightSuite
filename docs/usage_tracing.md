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
2. Select a fiber with `1`–`9` and click points around the **fiber edge** on each slice where it is visible.
3. Press `S` to save, then `F` to fit. Each fiber is fit as a circle/cylinder in CCF space and its sub-lens atlas regions are displayed.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** or **scroll** | Navigate coronal slices |
| **1 – 9** | Select the active fiber |
| **Click** | Add an edge point to the active fiber |
| **Backspace** | Delete the last point of the active fiber |
| **c** | Clear the active fiber's points on the current slice |
| **Enter** | Jump to a specific slice number |
| **s** | Save points to `grin_fiber_points.mat` |
| **f** | Fit all fibers and show atlas (and intensity) figures |

### Outputs

| File | Contents |
| :--- | :--- |
| `grin_fiber_points.mat` | Raw clicked edge points per fiber |
| `grin_fiber<N>_atlas.mat` | Per-fiber fit: `center_vox`, `normal_vox`, `radius_vox`, `atlas_pts`, `depths_um`, `slices_av`, `rvec_arr` |

Regenerate the figures from disk with `plotGRINTracingResults(savepath)`. The companion intensity figure is produced automatically when a `volume_registered/` folder is present.

### Tips

* Use **at least 4 edge points** per fiber, spread around the circumference and across several slices, so the circle fit is well constrained.
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
