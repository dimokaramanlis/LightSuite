# Slice Analysis

The Slice Module is optimized for conventional wide-field microscope data where you have individual coronal rather than a continuous 3D volume.

## Workflow Overview

The analysis is driven by the main script: `ls_analyze_slice_volume`.

The pipeline consists of the following stages:
1. **Volume Generation:** Load CZI files and create a centered 3D volume.
2. **Slice Reordering:** Manually reorder and flip slices as needed.
3. **Volume Alignment:** Initial rigid alignment to the Allen Atlas.
4. **Cutting Angle Determination:** Visually match the atlas cutting plane to your tissue.
5. **Control Point Matching:** Manually refine registration by placing landmarks.
6. **Registration Refinement:** Elastix-based affine and B-spline deformation per slice.
7. **Registered Volume Generation:** Apply transforms to all channels.
8. **Cell Detection:** Detect cells in the registered volume.
9. **Cell Registration to Atlas:** Transform cell coordinates to atlas space.

---

## Configuration Parameters

Parameters are stored in a `sliceinfo` struct defined at the top of `ls_analyze_slice_volume.m`.

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `sliceinfo.mousename` | Mouse identifier string | `'M001'` |
| `sliceinfo.channames` | Channel names | `{'DAPI', 'Cy3'}` |
| `sliceinfo.pxsize` | Original pixel size (µm) | `0.65` |
| `sliceinfo.px_process` | Processing resolution (µm) | `1.25` |
| `sliceinfo.px_register` | Registration resolution (µm) | `20` |
| `sliceinfo.px_atlas` | Atlas resolution (µm, typically 10) | `10` |
| `sliceinfo.slicethickness` | Physical spacing between slices (µm) | `40` |
| `sliceinfo.Nslices` | Total number of slices | `50` |
| `sliceinfo.celldiam` | Expected cell diameter (µm) | `14` |
| `sliceinfo.thresuse` | SNR thresholds `[detection, expansion]` | `[0.75, 0.4]` |
| `sliceinfo.debug` | Enable detection debug plots | `false` |
| `sliceinfo.use_gpu` | Use GPU for processing | `false` |
| `sliceinfo.atlasaplims` | Atlas AP axis limits `[min, max]` | `[200, 400]` |

---

## 1. Volume Generation

`generateSliceVolume(sliceinfo, regchan)` reads your CZI files and produces a centered 3D volume:

* Reads multi-scene CZI files using BioformatsImage.
* Applies a median filter to remove salt-and-pepper noise.
* Resamples to the processing resolution (`px_process`).
* Crops or pads each slice to a uniform size based on detected brain boundaries.
* Optionally denoises the DAPI channel.
* Saves a downsampled RGB preview (`volume_for_ordering.tiff`) for use in the next step.

---

## 2. Slice Reordering

Because slices may be loaded out of order or incorrectly oriented, `SliceOrderEditor` provides a GUI to review and correct the stack before proceeding.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** | Navigate to previous / next slice |
| **f** | Flip current slice horizontally |
| **o** | Toggle slice as excluded (removes from final stack) |
| **Enter** | Move current slice to a new position (prompts for index) |
| **s** | Save ordering decisions |
| **Escape** | Close the GUI |

**Output:** `volume_for_ordering_processing_decisions.txt` — a table recording the new order, flip state, and exclusion flags for each slice.

---

## 3. Volume Alignment

`alignSliceVolume(slicevol, sliceinfo)` performs an initial rigid alignment of the slice stack to the Allen Brain Atlas:

* Extracts a point cloud from each slice using difference-from-background (DFF) images.
* Loads the Allen template and annotations at the registration resolution.
* Runs a rigid → affine registration (`alignAtlasToSample`, `refineSampleFromAtlas`).

**Output:** `regopts.mat` — registration parameters including the rigid transform matrix (`tformrigid_allen_to_samp_20um`) and axis permutation.

---

## 4. Cutting Angle Determination

`determineCuttingAngleGUI(opts)` lets you visually match the 3D atlas cutting plane to the angle at which your tissue was sectioned.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** | Navigate slices |
| **Shift + ← / → / ↑ / ↓** | Rotate the 3D atlas view (0.3° per press) |
| **Scroll Wheel** | Move the atlas slice plane perpendicular to the view |
| **Spacebar** | Toggle atlas boundary overlay |
| **1 / 2 / 3 / 0** | Toggle individual channels (0 = all) |
| **Return** | Save the atlas position for the current slice |
| **c** | Clear the saved position for the current slice |

**Output:** `cutting_angle_data.mat` — camera direction vectors and atlas-space points defining the cutting plane for each slice.

---

## 5. Control Point Matching

`matchControlPointsInSlices(opts)` opens a GUI to manually refine the 2D registration of each slice by placing corresponding anatomical landmarks.

### Interface
* **Left panel:** Your histology slice.
* **Right panel:** The corresponding atlas plane.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **Click (Left panel)** | Place a control point on the histology image |
| **Click (Right panel)** | Place the corresponding control point on the atlas |
| **← / →** | Navigate to previous / next slice |
| **Spacebar** | Toggle atlas overlay |
| **g** | Toggle grid overlay |
| **c** | Clear control points for the current slice |
| **1 / 2 / 3 / 0** | Toggle channels |
| **Enter** | Jump to a specific slice number |
| **s** | Save control points |

**Tips:**
* At least 3 matched point pairs are required per slice for an affine fit.
* Spread points across the full extent of the slice — don't cluster them in one region.
* The MSE displayed updates live; lower is better.

**Output:** `atlas2histology_tform.mat` — per-slice affine transforms and control point arrays.

---

## 6. Registration Refinement

`registerSlicesToAtlas(opts)` applies a two-stage registration per slice:

1. **Affine stage:** Uses your control points (if ≥ 5 are available) or falls back to image-based affine fitting.
2. **B-spline stage:** Runs an Elastix deformable registration to capture local deformations not covered by the affine.

Both forward (atlas → sample) and reverse (sample → atlas) transforms are saved, as the reverse is needed for mapping cell coordinates.

**Outputs:**
* `transform_params.mat` — complete per-slice registration parameters.
* `elastix_forward/NNN_slice_bspline_atlas_to_samp_20um.txt`
* `elastix_reverse/NNN_slice_bspline_samp_to_atlas_20um.txt`

---

## 7. Registered Volume Generation

`generateRegisteredSliceVolume(sliceinfo, transformparams)` maps all channels into atlas space:

* For each slice, applies the inverse B-spline deformation followed by the inverse affine transform.
* Places each transformed slice at its corresponding atlas Z position.
* Fills gaps between slices via nearest-neighbor interpolation.
* Applies the global rigid transform to place the volume in Allen atlas coordinates.

**Output:** `volume_registered` — the full 3D volume in atlas space for all channels.

---

## 8. Cell Detection

`extractCellsFromSliceVolume(opts, ichan)` detects cells slice by slice using a 2D SNR-based approach.

### Algorithm
1. **Background subtraction:** Computes a local background using a median filter and generates a difference-from-background (DFF) image.
2. **Bandpass filtering:** Enhances objects matching the expected cell size while suppressing noise and background.
3. **Local maxima detection:** Identifies candidate cells as local maxima above `thresuse(1)`.
4. **Morphological filtering:** Removes candidates that are too elongated (circularity < 0.7), too small, or too dim (below `thresuse(2)`).

### Key Parameters
* **`celldiam`**: Expected cell diameter in µm. This directly controls the bandpass filter kernel size — it is the most important detection parameter.
* **`thresuse`**: Two-element SNR threshold vector.
  * `thresuse(1)` — initial detection threshold. Lower values detect dimmer cells but increase false positives.
  * `thresuse(2)` — expansion threshold used for cell boundary growth and minimum intensity filtering.

**Output:** `cell_locations_sample.mat` — an N×5 array with columns `[x, y, slice_index, diameter, mean_intensity]`.

If `debug = true`, PNG overlays are saved to `cell_detections/` showing the DFF image and detected cell masks per slice.

---

## 9. Cell Registration to Atlas

`slicePointsToAtlas(inputpts, transformparams)` transforms detected cell coordinates from slice space into Allen Atlas space via a chain of transforms:

1. B-spline deformation field (elastix transformix)
2. Inverse affine per slice
3. Global inverse rigid transform

**Output:** `cell_locations_atlas.mat` — an N×3+ array with columns `[ML, AP, DV]` in atlas voxels, plus the original diameter and intensity properties.

---

## Output File Summary

```
<save_folder>/
├── volume_centered                    # Full-resolution centered slice volume
├── volume_for_ordering.tiff           # RGB preview for the ordering GUI
├── volume_for_ordering_processing_decisions.txt
├── volume_ordered.tiff                # Reordered/flipped volume
├── sample_register_20um.tif           # Volume at registration resolution
├── volume_for_inspection.tiff         # Downsampled RGB for inspection
├── regopts.mat                        # Initial alignment parameters
├── sliceinfo.mat                      # Slice metadata
├── cutting_angle_data.mat             # Cutting plane angles/vectors
├── atlas2histology_tform.mat          # Control point registration data
├── transform_params.mat               # Final per-slice registration parameters
├── elastix_forward/                   # Forward B-spline transforms
│   └── NNN_slice_bspline_atlas_to_samp_20um.txt
├── elastix_reverse/                   # Reverse B-spline transforms (for cell mapping)
│   └── NNN_slice_bspline_samp_to_atlas_20um.txt
├── volume_registered                  # Final atlas-space volume (all channels)
├── cell_locations_sample.mat          # Detected cells in sample space
├── cell_locations_atlas.mat           # Detected cells in atlas space
└── cell_detections/                   # Debug visualizations (if debug=true)
    └── slice_NNN_detections.png
```
