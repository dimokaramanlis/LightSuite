# Slice Analysis

The slice module turns serial wide-field sections into a 3D volume in Allen CCFv3 space. You order the sections, align them rigidly to the atlas, optionally correct the cutting angle, and register each section non-rigidly with the same landmark-guided objective used for volumes. Gaps between sections are filled by nearest-neighbour interpolation.

> Shared concepts (registration, control points, cell detection) are explained in [How it works](how_it_works.md). This page covers the slice-specific steps.

Run the sections of `demos/ls_analyze_slice_volume.m` in order. Sections marked *(manual)* open a GUI.

| Step | Function | Type |
| :--- | :--- | :--- |
| 1. Build the slice volume | `generateSliceVolume` | auto |
| 2. Order, flip, exclude and crop slices | `SliceOrderEditor`, `generateReordedVolume` | manual |
| 3. Rigid alignment to the atlas | `alignSliceVolume` | auto |
| 4. Cutting angle (optional) | `determineCuttingAngleGUI` | manual |
| 5. Control points | `matchControlPointsInSlices` | manual |
| 6. Per-section registration | `registerSlicesToAtlas` | auto |
| 7. Registered volume | `generateRegisteredSliceVolume` | auto |
| 8. Cell detection | `extractCellsFromSliceVolume` | auto |
| 9. Cells to atlas | `slicePointsToAtlas` | auto |
| 10. Region tables | `generateSliceIntensitiesCSV`, `sliceCellCountsToCSV` | auto |

---

## Input and settings

Set `datafolderpath` and `mousename` at the top of the demo. `datafolderpath` can hold one subfolder per mouse (matched by `mousename`) or be the image folder itself. Accepted input:

* **CZI** — each file can hold several scenes (one per slice); pixel size and channel names are read from the metadata.
* **TIF/TIFF** — one file per slice, channels as pages. Set `sliceinfo.pxsize` and `sliceinfo.channames` in the demo.

Settings are read from `local_settings.txt` in the sample folder. If it is missing, it is created with the defaults below; edit it and re-run the section. The file is read again before alignment (step 3), so later edits are picked up.

| Setting | Default | Meaning |
| :--- | :--- | :--- |
| `slicethickness` | `150` | Spacing between consecutive sections (µm) |
| `px_process` | `5` | Processing resolution (µm/px) |
| `px_register` | `20` | Registration resolution (µm/px) |
| `px_atlas` | `10` | Atlas resolution (µm/px) |
| `atlasaplims` | `[180 1079]` | AP range of the atlas used (10 µm voxels) |
| `regchan` | `'dapi'` | Channel used for registration (name substring) |
| `denoisedapi` | `false` | Local-std filter on the DAPI channel |
| `medianfiltreg` | `false` | Median filter on the registration volume |
| `use_gpu` | `false` | Use the GPU where supported |

All output goes to a `lightsuite/` folder next to the images.

---

## 1. Build the slice volume

`generateSliceVolume` median-filters each image (salt-and-pepper noise), resamples it to `px_process`, pads all slices to a common size, and puts the registration channel first. It also detects the brain outline in each slice and stores a crop suggestion, which is applied in step 3.

**Outputs:** `volume_centered/`, `volume_for_ordering.tiff` (RGB preview at `px_register`), `volume_for_ordering_processing_decisions.txt` (crop suggestions), `sliceinfo.mat`.

## 2. Order slices

`SliceOrderEditor` shows the preview stack so you can fix order, orientation and bad slices. `generateReordedVolume` then writes the result to `volume_ordered.tiff` for a quick check.

| Key | Action |
| :--- | :--- |
| **← / →** | Previous / next slice |
| **Enter** | Move the current slice to a new position |
| **f** | Flip the slice left-right |
| **o** | Exclude / restore the slice |
| **c** | Edit the crop rectangle (double-click to confirm, Esc to cancel) |
| **m** | Show / hide the montage (click a tile to jump to it, drag it onto another tile to move it) |
| **1 / 2 / 3** | Show / hide channel |
| **d** | Show / hide the display panel (channel levels; display only) |
| **s** | Save decisions |
| **Esc** | Close (asks to save) |

**Outputs:** `volume_for_ordering_processing_decisions.txt` (order, flip/exclude state, crop per slice); `slice_order_montage.png` if the montage is open when saving.

## 3. Rigid alignment

`alignSliceVolume` applies the ordering and crops, then aligns the stack to the atlas:

* A point cloud of bright pixels is taken from each section and stacked along AP using `slicethickness`.
* Three passes alternate a global 3D rigid fit to the atlas (Coherent Point Drift) with per-section 2D rigid fits (ICP) against the atlas plane at each section.

**Check the diagnostic plots and the inspection volume before going on** — the control point GUI starts from this alignment.

**Outputs:** `regopts.mat` (rigid transform `tformrigid_allen_to_samp_20um`), `volume_aligned/chanNN_<name>.tiff` (all channels, `px_process`), `sample_register_20um.tif`, `volume_for_inspection.tiff`, `<mousename>_dim{1,2,3}_initial_registration.png`.

## 4. Cutting angle (optional)

Use this if sections were not cut orthogonal to the AP axis. `determineCuttingAngleGUI` lets you rotate the atlas until its plane matches a section and save that direction. Saved directions are averaged into one cutting normal, which corrects the rotation of the rigid transform in steps 5 and 6 — a handful of sections is enough.

| Key | Action |
| :--- | :--- |
| **← / →** | Previous / next section |
| **Shift + arrows** | Rotate the atlas (0.3° per press) |
| **Scroll** | Move the atlas plane along its normal |
| **Space** | Toggle atlas boundaries |
| **1 / 2 / 3, 0** | Show one channel / all channels |
| **Enter** | Store the atlas plane for this section |
| **c** | Clear the stored plane for this section |
| **s** | Save (also offered on close) |

**Output:** `cutting_angle_data.mat`. Delete it to fall back to the automatic angle.

## 5. Control points

`matchControlPointsInSlices` shows each section next to its predicted atlas plane. Click matching landmarks on both; the title shows the number of pairs and the error of a live affine fit (3 pairs minimum). The atlas AP position of sections you have not annotated is predicted from those you have, so annotating a few sections spread along the stack also fixes spacing errors.

| Key | Action |
| :--- | :--- |
| **Click** | Add a control point |
| **Backspace** | Delete the last point |
| **← / →** | Previous / next section |
| **Scroll** (Shift: ×10) | Move the atlas plane along AP |
| **Space** | Toggle atlas boundaries |
| **g** | Toggle grid |
| **c** | Clear this section's points |
| **1 / 2 / 3, 0** | Show one channel / all channels |
| **s** | Save (also offered on close) |

See [How it works → control points](how_it_works.md#refining-registration-with-control-points-active-learning) for where to place landmarks.

**Output:** `atlas2histology_tform.mat`.

## 6. Per-section registration

`registerSlicesToAtlas(opts)` registers each section in two stages:

1. **Affine** — fitted to your control points if a section has at least 5, otherwise to image-derived point clouds.
2. **B-spline** (elastix) — Mattes mutual information plus a control-point distance term weighted by `opts.cpwt` (demo: 0.2; increase for damaged sections).

The B-spline grid spacing adapts per section: half the square root of the brain area without points, refining toward `opts.bspline_spatial_scale` (default 0.96 mm) as the points cover the section more evenly. Clustered points barely refine it.

Sections without control points take their AP position from the annotated ones (interpolated when more than 3 sections are annotated).

**Outputs:** `transform_params.mat`; `elastix_forward/NNN_slice_bspline_atlas_to_samp_20um.txt` (used to map points) with `NNN_slice_registration_comparison.png`; `elastix_reverse/NNN_slice_bspline_samp_to_atlas_20um.txt` (used to warp images).

## 7. Registered volume

`generateRegisteredSliceVolume` warps every channel of each section into its atlas plane (inverse B-spline, then inverse affine), fills the gaps between sections by nearest-neighbour interpolation along AP, and places the stack in the atlas with the inverse rigid transform.

**Output:** `volume_registered/` — all channels at 10 µm in atlas space.

## 8. Cell detection

`extractCellsFromSliceVolume(sliceinfo, ichans)` runs a 2D version of the [SBR detector](how_it_works.md#cell-detection-and-artifact-classification) on each section of `volume_aligned`, for one or more channels. The background is a local rank filter (2.5 cell diameters wide); candidates with circularity below 0.7 are discarded. Set these in the demo:

| Parameter | Demo | Meaning |
| :--- | :--- | :--- |
| `celldiam` | `14` | Expected cell diameter (µm); sets the filter scale — the most important parameter |
| `thresuse` | `[0.75 0.4]` | SBR thresholds: detection, then expansion / minimum mean intensity |
| `debug` | `true` | Save detection overlays to `cell_detections/` (slower) |

**Output:** `chanNN_cell_locations_sample.mat` — `cell_locations`, N×5 `[x, y, section, diameter, intensity]`, in `volume_aligned` pixels. View with `visualizeCellDetections(procpath, Space = 'sample')`.

## 9. Cells to atlas

`slicePointsToAtlas(cell_locations, transformparams)` maps cells through the forward B-spline, the inverse affine and the inverse rigid transform — the exact reverse of step 7, so a cell lands on the voxel that shows it in `volume_registered`.

**Output:** `chanNN_cell_locations_atlas.mat` — `atlasptcoords`, `[DV, AP, ML]` in 10 µm atlas voxels, then diameter and intensity. View with `visualizeCellDetections(procpath, Space = 'atlas')`.

## 10. Region tables

* `generateSliceIntensitiesCSV(sliceinfo, 'writetocsv', true)` — median intensity per region and hemisphere, per channel.
* `sliceCellCountsToCSV(procpath, 'writetocsv', true)` — cell counts per region and hemisphere, per channel.

**Outputs** (in `volume_registered/`): `chanNN_intensities.mat/.csv`, `chanNN_cellcounts.mat/.csv`.

---

## Output files

```
<sample folder>/
├── local_settings.txt
└── lightsuite/
    ├── sliceinfo.mat
    ├── volume_centered/                        # step 1, unordered
    ├── volume_for_ordering.tiff
    ├── volume_for_ordering_processing_decisions.txt
    ├── slice_order_montage.png
    ├── volume_ordered.tiff
    ├── volume_aligned/chanNN_<name>.tiff       # ordered, cropped, rigidly aligned
    ├── volume_for_inspection.tiff
    ├── sample_register_20um.tif
    ├── <mousename>_dim{1,2,3}_initial_registration.png
    ├── regopts.mat
    ├── cutting_angle_data.mat                  # optional
    ├── atlas2histology_tform.mat
    ├── transform_params.mat
    ├── elastix_forward/                        # point mapping + QC images
    ├── elastix_reverse/                        # image warping
    ├── volume_registered/                      # atlas-space volume + region tables
    ├── chanNN_cell_locations_sample.mat
    ├── chanNN_cell_locations_atlas.mat
    └── cell_detections/                        # if debug = true
```
