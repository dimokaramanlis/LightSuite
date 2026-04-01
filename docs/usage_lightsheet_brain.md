# Lightsheet Whole-Brain Analysis

This module processes whole-brain datasets acquired via lightsheet microscopy. It handles large volumetric data (100 GB+), performing preprocessing, cell detection, and registration to the Allen Brain Atlas (CCF v3).

## Before You Start

You will need:
* A stitched TIFF dataset (output from Terastitcher or BigStitcher).
* A fast SSD with at least 500 GB of free space for temporary processing files (`opts.fproc`).
* [Elastix](https://elastix.lumc.nl/) installed and on your system path (used for deformable registration).

---

## Workflow Overview

The analysis is driven by the main script: `ls_analyze_lightsheet_volume`.

The pipeline consists of six stages, some automated and some requiring your input:

| Step | Stage | Type |
| :---: | :--- | :--- |
| 1 | Preprocessing | Automated |
| 2 | Initial Registration | Automated |
| 3 | Manual Alignment (Control Point GUI) | **Manual** |
| 4 | Deformable Registration | Automated |
| 5 | Volume Registration to Atlas | Automated |
| 6 | Cell Detection & Atlas Mapping | Automated |

---

## Configuration Parameters

Open `ls_analyze_lightsheet_volume.m` and fill in the `opts` struct at the top of the script before running anything.

### Paths

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `opts.mousename` | A short name for this sample | `'M001'` |
| `opts.datafolder` | Path to your stitched TIFF files | `'D:\Data\M001'` |
| `opts.fproc` | Path to fast SSD for temporary binary files (≥500 GB) | `'E:\proc\M001'` |
| `opts.savepath` | Where results will be saved | `fullfile(opts.datafolder, 'lightsuite')` |

### Data Format

| Parameter | Description | Options |
| :--- | :--- | :--- |
| `opts.tifftype` | How your TIFF files are organized | `'channelperfile'` (BigStitcher) or `'planeperfile'` (Terastitcher) |
| `opts.pxsize` | Voxel size in microns `[x y z]` | `[20 20 20]` |

### Cell Detection

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `opts.channelforcells` | Which channel contains your labeled cells | `3` |
| `opts.celldiam` | Approximate cell diameter in microns | `25` |
| `opts.thres_cell_detect` | SNR thresholds `[detection, expansion]` | `[0.5, 0.4]` |
| `opts.savecellimages` | Save 2D projections of each detected cell | `false` |

> **Tip:** `celldiam` is the most important detection parameter. Measure a few representative cells in ImageJ and use that value. If you get too many false positives, raise `thres_cell_detect(1)`.

### Registration

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `opts.channelforregister` | Which channel to use for atlas registration (use an autofluorescence or structural channel) | `2` |
| `opts.bspline_spatial_scale` | B-spline deformation scale in mm. Smaller = finer local warping | `0.64` |
| `opts.augmentpoints` | Automatically add detected landmarks to supplement manual control points | `false` |
| `opts.weight_usr_pts` | How much weight to give your manual control points relative to image data | `0.2` |

### Output Options

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `opts.debug` | Save diagnostic images for cell detection | `false` |
| `opts.writetocsv` | Export regional intensity and cell count tables to CSV | `false` |

---

## 1. Preprocessing

**Function:** `preprocessLightSheetVolume(opts)`

This stage reads your raw TIFF files and prepares two versions of the data:

* **Full-resolution binary** — written to `opts.fproc` for cell detection. The raw data is filtered with a 3×3 median filter to remove salt-and-pepper noise, then saved as a binary file for fast random access during detection.
* **Downsampled registration volume** — a lower-resolution TIFF (`chan_X_sample_register_20um.tif`) used for all registration steps. This is also median-filtered and rescaled to 16-bit.

**Outputs:**
* `chan_X_sample_register_20um.tif` — one per channel, in `opts.savepath`
* `chan_X_binary_MOUSENAME.dat` — full-resolution binary, in `opts.fproc`
* `regopts.mat` — saves the current `opts` struct for later stages

---

## 2. Initial Registration

**Function:** `initializeRegistration(opts.savepath)`

This stage automatically finds a coarse alignment between your sample and the Allen Atlas before you refine it manually:

1. Loads the downsampled registration volume and the Allen CCF template.
2. Prompts you (if needed) to set the **brain orientation** — which axis is anterior-posterior, which is dorsal-ventral, and whether any axis needs to be flipped. This is saved to `brain_orientation.txt` so you only do it once.
3. Extracts point clouds from both volumes based on tissue edges.
4. Runs a similarity transform (rotation + scale) to coarsely align the sample to the atlas.
5. Auto-detects candidate anatomical landmarks as a starting point for the next step.

**Outputs:**
* `brain_orientation.txt`
* `regopts.mat` (updated with initial transform and auto-detected landmarks)
* `dim1/2/3_initial_registration.png` — diagnostic images showing the initial fit

---

## 3. Manual Alignment (Control Point GUI)

**Function:** `matchControlPoints_unified(opts)`

This is the main step requiring your attention. The GUI shows your sample alongside the corresponding Allen Atlas slice so you can click matching anatomical landmarks in both views.

### Interface
* **Left panel:** Your histology/sample slice.
* **Right panel:** The corresponding Allen Atlas plane.
* **Top banner:** Current fit quality (MSE) and total point count.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **Click (Left panel)** | Place a control point on the sample |
| **Click (Right panel)** | Place the corresponding point on the atlas |
| **← / →** | Navigate to the previous / next slice |
| **Scroll Wheel** | Move the atlas slice plane independently |
| **Backspace** | Delete the last point added |
| **c** | Clear all points on the current slice |
| **Spacebar** | Toggle the red atlas boundary overlay |
| **Enter** | Jump to a specific slice number |
| **s** | Save and exit |

### Tips for Good Results
* **Use at least 16 matched pairs** spread across the full AP extent of the brain. The fit only activates once 16 points are placed. 100 points or more are recommended.
* **Slices are shown in randomized order.** This is intentional — it encourages you to spread points evenly rather than clustering them in one region.
* **Watch the MSE** in the top banner. It updates live as you add points. A lower MSE means a tighter correspondence and less locally deformed tissue.
* **Good landmarks** include the boundaries of ventricles, major fiber tracts (corpus callosum, anterior commissure), and distinct nuclei outlines. Avoid regions with tissue damage.
* If a region is damaged in your sample, estimate where the structure *would* be and place a point there anyway. This prevents the registration from bending the atlas into the damaged area.

**Output:** `atlas2histology_tform.mat` — affine transform and all control point arrays.

---

## 4. Deformable Registration

**Function:** `multiobjRegistration(opts)`

This stage refines the alignment further using your control points and a deformable B-spline registration:

1. **Affine step:** Fits a global affine transform using your control points combined with auto-detected landmarks.
2. **B-spline step:** Runs Elastix to compute a smooth, local deformation field that captures non-linear tissue distortion.
3. **Inverse transform:** Computes the reverse mapping (sample → atlas), needed for cell coordinate mapping.

**Outputs:**
* `transform_params.mat` — the complete registration parameters
* `elastix_forward/bspline_atlas_to_samp_20um.txt`
* `elastix_reverse/bspline_samp_to_atlas_20um.txt`
* `*_affine_registration.png` / `*_bspline_registration.png` — diagnostic images

---

## 5. Volume Registration to Atlas

**Function:** `generateRegisteredBrainVolumes(opts.savepath)`

Applies the registration to all channels and computes regional statistics:

* Warps each channel's registration volume into 10 µm Allen Atlas space using the B-spline transform.
* For each brain region and hemisphere, computes the **median signal intensity** and **volume in mm³**.

**Outputs** (in `volume_registered/`):
* `chan0X_intensities.mat` — median intensity and volume per region, per hemisphere
* `chan0X_intensities.csv` — same data as a table (if `writetocsv = true`), with columns: `name`, `structure`, `division`, `parcellation_index`, `RightIntensity`, `LeftIntensity`, `RightVolume[mm3]`, `LeftVolume[mm3]`

---

## 6. Cell Detection & Atlas Mapping

Cell detection runs automatically during preprocessing (Step 1) on the channel specified by `opts.channelforcells`. The detected cell coordinates are then mapped to atlas space.

### Detection Algorithm

Detection is performed in 3D batches on the full-resolution data:

1. **3D Bandpass Filtering:** Enhances signal at the expected cell size while suppressing noise and background.
2. **Background Estimation:** Estimates the local noise floor from a statistical sample of voxels.
3. **Local Maxima Detection:** Finds candidate cells as local maxima above `thres_cell_detect(1)`.
4. **Morphological Filtering:** Removes candidates that are too elongated (axis ratio > 2.5), too small, or too dim.

### Key Parameters
* **`celldiam`** — controls the bandpass filter. Set this to your actual cell diameter.
* **`thres_cell_detect(1)`** — primary detection SNR threshold. Higher = fewer, more confident detections.
* **`thres_cell_detect(2)`** — secondary threshold used during cell boundary expansion and minimum intensity filtering.

### Atlas Mapping

**Function:** `transformPointsToAtlas(opts.savepath)`

Detected cells are transformed to atlas space through the same chain of transforms used for the volume: similarity → affine → B-spline.

**Outputs** (in `volume_registered/`):
* `chan_X_cell_locations_atlas.mat` — N×6 array with columns `[x, y, z, intensity, diameter, ellipticity]` in atlas voxel coordinates
* `cell_counts_by_region.csv` — cell counts and median intensities per region and hemisphere (if `writetocsv = true`)

---

## Output File Summary

```
<savepath>/
├── regopts.mat                              # Configuration (updated at each stage)
├── brain_orientation.txt                    # Axis permutation (set once)
├── atlas2histology_tform.mat                # Manual control point alignment
├── transform_params.mat                     # Final registration parameters
├── chan_1_sample_register_20um.tif          # Downsampled volumes (per channel))
├── chan_1_cell_locations_sample.mat         # Detected cells in sample space
├── chan_1_cell_detections/                  # Debug images (if debug=true)
├── dim1/2/3_initial_registration.png        # Initial alignment check
├── *_affine_registration.png                # After affine step
├── *_bspline_registration.png               # After B-spline step
└── volume_registered/
    ├── chan01_intensities.mat               # Regional intensity & volume stats
    ├── chan01_intensities.csv               # (if writetocsv=true)
    ├── chan_1_cell_locations_atlas.mat      # Cells in atlas coordinates
    └── cell_counts_by_region.csv           # Regional cell counts (if writetocsv=true)
```
