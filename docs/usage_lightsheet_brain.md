# Whole-Brain Volume Analysis

This module registers whole-brain volumes to the Allen CCFv3 and detects labelled cells at native resolution. It handles 100 GB+ datasets on a single workstation, yielding ~4.5 × 10⁵ cells per brain in 3–4 h.

**Any volumetric modality works.** The pipeline sees a stitched 3D stack, not a microscope: light-sheet, fMOST, serial two-photon tomography and comparable volumetric datasets all go through unchanged. Set `opts.pxsize` to your voxel size and `opts.celldiam` to your cell size; nothing else is modality-specific. Light-sheet is simply what the reference figures were acquired with.

> New to LightSuite? Read [How it works](how_it_works.md) first — it explains the registration stages, control-point active learning, and SBR cell detection that this page builds on. Here we cover the brain-specific steps and parameters.

## Before You Start

You will need:
* A stitched TIFF dataset (output from Terastitcher or BigStitcher).
* A fast SSD with at least 500 GB of free space for temporary processing files (`opts.fproc`).
* [Elastix](https://elastix.lumc.nl/) installed and on your system path (used for deformable registration).

---

## Workflow Overview

The analysis is driven by the main script: `ls_analyze_lightsheet_volume`.

The pipeline consists of seven stages, some automated and some requiring your input:

| Step | Stage | Type |
| :---: | :--- | :--- |
| 1 | Preprocessing | Automated |
| 2 | Initial Registration | Automated |
| 3 | Manual Alignment (Control Point GUI) | **Manual** |
| 4 | Deformable Registration | Automated |
| 5 | Volume Registration to Atlas | Automated |
| 6 | Cell Detection & Atlas Mapping | Automated |
| 7 | CNN Cell Classification | Optional, one-off **manual** labelling |

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
| `opts.pxsize` | Voxel size in microns `[x y z]` | `[6.55 6.55 5]` |

### Cell Detection

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `opts.channelforcells` | Which channel contains your labeled cells | `3` |
| `opts.celldiam` | Approximate cell diameter in microns | `14` |
| `opts.thres_cell_detect` | SBR thresholds `[primary, secondary]` | `[0.5, 0.4]` |
| `opts.savecellimages` | Save 2D projections of each detected cell | `false` |

> **Tip:** `celldiam` is the most important detection parameter. Measure a few representative cells in ImageJ and use that value. If you get too many false positives, raise `thres_cell_detect(1)`.

> **Set `savecellimages = true` if you plan to use the [CNN classifier](#7-cnn-cell-classification).** The three-view images it needs are written only during detection; without them the classifier cannot be applied afterwards and the volume has to be re-processed.

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

See [How it works → control points](how_it_works.md#refining-registration-with-control-points-active-learning) for landmark strategy. In short: the fit activates at 16 pairs and keeps improving to about 100, which takes 10–20 minutes; watch the live MSE; favour ventricles and fiber tracts; and in damaged regions place a point where the structure *should* be. Views are shown in randomized order to encourage even coverage.

**Output:** `atlas2histology_tform.mat` — affine transform and all control point arrays.

---

## 4. Deformable Registration

**Function:** `multiobjRegistration(opts)`

This stage refines the alignment further using your control points and a deformable B-spline registration:

1. **Affine step:** Fits a global affine transform by least squares from your control points plus the auto-detected correspondences. Automatic points within 1 mm of a user point are dropped so they cannot dilute yours.
2. **B-spline step:** Elastix computes a free-form deformation on a 0.64 mm control-point grid, optimizing Advanced Mattes mutual information (48 histogram bins, weight 1.0) *plus* a corresponding-points distance term (weight `opts.weight_usr_pts`, default 0.2), over four resolution levels with adaptive stochastic gradient descent. The coarse grid is what prevents overfitting — it matches the several-hundred-µm scale at which fixation and clearing actually deform tissue.
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

### Choosing the per-area statistic

The median is the default because it ignores the bright outliers (vessels, debris, labelled somata) that would drag a mean around. Pass `'areafun'` to summarize each area some other way:

```matlab
generateRegisteredBrainVolumes(opts.savepath, 'areafun', @mean);          % mean intensity
generateRegisteredBrainVolumes(opts.savepath, 'areafun', @std);           % spread within each area
generateRegisteredBrainVolumes(opts.savepath, 'areafun', @(x) quantile(x, 0.9));
```

Any function handle that takes a vector of voxel values and returns one number works; it is checked on a test vector before the loop starts, so a wrong handle fails immediately rather than after the warping. The same function is applied to the out-of-brain background level, so the two stay comparable. You can also set it once as `opts.areafun` in `regopts.mat`.

The per-channel `.mat` keeps the variable name `medianoverareas` whatever statistic you choose — existing readers such as `loadMouseBackgroundSignal` keep working — and records the function used as `areafun` / `areafunname` alongside it.

**Outputs** (in `volume_registered/`):
* `chan0X_intensities.mat` — median intensity and volume per region, per hemisphere
* `chan0X_intensities.csv` — same data as a table (if `writetocsv = true`), with columns: `name`, `structure`, `division`, `parcellation_index`, `RightIntensity`, `LeftIntensity`, `RightVolume[mm3]`, `LeftVolume[mm3]`

---

## 6. Cell Detection & Atlas Mapping

Cell detection runs automatically during preprocessing (Step 1) on the channel specified by `opts.channelforcells`. The detected cell coordinates are then mapped to atlas space.

### Detection Algorithm

Detection runs on the full-resolution channel in overlapping GPU 3D batches (1800 × 1800 × 32 voxels): band-pass filter → signal-to-background ratio → local maxima → morphological filter → CNN artifact classifier. See [How it works → cell detection](how_it_works.md#cell-detection-and-artifact-classification) for the full description.

### Key Parameters
* **`celldiam`** — controls the band-pass filter. Set this to your actual cell diameter.
* **`thres_cell_detect(1)`** — primary SBR threshold for accepting a local maximum. Higher = fewer, more confident detections.
* **`thres_cell_detect(2)`** — secondary threshold used during cell boundary expansion and minimum intensity filtering.

### Atlas Mapping

**Function:** `transformPointsToAtlas(opts.savepath)`

Detected cells are transformed to atlas space through the same chain of transforms used for the volume: similarity → affine → B-spline.

**Outputs** (in `volume_registered/`):
* `chan_X_cell_locations_atlas.mat` — N×6 array with columns `[x, y, z, intensity, diameter, elongation]` in atlas voxel coordinates
* `cell_counts_by_region.csv` — cell counts and median intensities per region and hemisphere (if `writetocsv = true`)

### Quantification conventions

Quantification happens in atlas space, at the **substructure level** of the Allen parcellation, **separately per hemisphere** (the midplane splits the cells):

* Detections outside the brain annotation or beyond the volume bounds are discarded.
* Regional **volume** is the annotation voxel count × voxel volume (10⁻⁶ mm³ at 10 µm); **density** is count ÷ volume.
* Regional **intensity** is the median voxel value in the region, reported **relative to the median out-of-brain intensity** of that hemisphere. Use `'areafun'` to summarize with something other than the median.

> **Soma diameter is a relative measure.** The `diameter` column comes from the SBR-dilated candidate region, so it overestimates true anatomy (~20 µm reported against 10–16 µm actual), partly through the dilation and partly through point-spread-function broadening in undeconvolved data. Laminar and cross-region *comparisons* are faithful; absolute values are not.

For cohort work, note that labelling efficiency varies several-fold between animals — normalize before pooling rather than averaging raw densities.

### Bringing in points from elsewhere

Cells counted outside LightSuite can go through the same transform. Three formats are accepted, and a single file can be passed directly with `'savepath'` pointing at the registration folder:

| Format | Contents | Cell images? |
| :--- | :--- | :--- |
| `.mat` | a `cell_locations` array `[N × M]`, `M ≥ 3`, columns `[x y z, descriptors…]` | yes, if saved |
| `.csv` | the same array as text — what `writematrix(cell_locations, …)` writes. A header line is allowed and skipped | no |
| `.xml` | an ImageJ / Fiji **Cell Counter** marker file; each `<Marker_Type>` block becomes one point set | no |

```matlab
% one ImageJ Cell Counter file, markers counted on channels 2 and 3
transformPointsToAtlas('D:\data\cells.xml', 'savepath', opts.savepath, ...
    'channel', [2 3], 'writetocsv', true);
```

In folder mode the search patterns are `*cell_locations_sample.mat`, `*cell_locations_sample.csv` and `*.xml`; an XML that turns out not to be a Cell Counter file is skipped with a warning rather than aborting the run.

**Every point set has to name its channel**, because the regional statistics are written per channel (`chanNN_cellcounts.mat`). The channel is resolved in this order:

1. the `'channel'` option — a scalar for all point sets, or one entry per point set;
2. a `chan_3_` / `chan03_` / `channel3_` tag in the file name;
3. for XML only, the `<Type>` of the marker block.

A file that names no channel at all falls back to its position in the list, as it always has, but **warns that it guessed** — rename it to `chan_<N>_...` or pass `'channel'` to make it definite. The mapping it settled on is printed for every point set either way, and two point sets claiming the same channel produce a warning, since their statistics would overwrite each other.

> **Marker types are not channels.** The Cell Counter plugin's `<Type>` is a counter category; nothing in the file records which channel it was counted on. Falling back to it is a convenience for the common case where they line up — pass `'channel'` whenever they do not.

> **XML coordinate convention.** The plugin stores `MarkerX`/`MarkerY` as 0-based canvas pixels and `MarkerZ` as the 1-based ImageJ slice number, while LightSuite's `cell_locations` are 1-based throughout. `[1 1 0]` is therefore added to every marker. Pass `'coordoffset', [0 0 0]` to take the file's numbers exactly as they are.

---

## 7. CNN Cell Classification

Detection is deliberately permissive — it is easier to reject a bad candidate than to recover a missed cell — so the raw output contains false positives: bubble edges, tissue–solution interfaces, bright neurite segments, tissue folds. A small convolutional network trained on your own labels removes them.

The network sees each candidate as **three maximum-intensity views** (`xy`, `xz`, `yz`) stacked as the channels of one 37×37×3 image. Classifying three projections rather than a 3-D volume keeps it small enough to run over a million candidates in minutes while still capturing the shape asymmetries that separate cells from artifacts.

### Step 1 — detect with cell images

The classifier needs the views, and they are written **only during detection**:

```matlab
opts.savecellimages = true;
```

`chan_X_cell_locations_sample.mat` then also holds `cell_images` (the three flattened views) and `imwindow` (their half-window). Without them the volume has to be re-processed.

### Step 2 — label a training set

**Function:** `CellLabelingTool`

Opens a file picker, loads a detection `.mat`, and shows the candidates one at a time in a random order.

| Key | Action |
| :--- | :--- |
| **C** | Label as **cell** (1) |
| **N** | Label as **noise** (0) |
| **U** | Unlabel (back to unsorted) |
| **← / →** | Back / next (skip) |
| **S** | Save progress |

Progress is saved to `<name>_labeled.mat` and the tool resumes on reopen, so labelling can be spread over sittings. A few hundred candidates per brain is usually enough; label from more than one brain if your staining varies. Only labelled candidates are saved, so skipping the ambiguous ones is fine.

### Step 3 — train

**Demo script:** `demos/trainClassificationNetwork.m`

Point `textfilewithpaths` at a text file listing your `*_labeled.mat` files, one per line, and set `netowrksavepath`. The script pools all of them, balances the two classes (`balanceChoices`), splits 80/20 into training and validation, and trains a small CNN — three conv/batchnorm/ReLU blocks with pooling, then a two-way fully connected layer — with rotation, reflection and translation augmentation (`augmentCellViews`) on the training half only.

It saves `<date>_CellClassifierNet.mat` containing `net` plus the training, validation and overall accuracies, and calls `visualizeNetworkClassification` so you can see what it got wrong. The classifier reported in the paper reaches **98.7%** validation accuracy and discards 10–20% of candidates (median 16%); the networks bundled here sit around 97%.

### Step 4 — apply during atlas mapping

**Function:** `transformPointsToAtlas(..., 'network', NET)`

```matlab
transformPointsToAtlas(opts.savepath, 'writetocsv', true, ...
    'network', 'D:\nets\20260602_CellClassifierNet.mat');
```

`NET` is a path to the saved `.mat` or the network object itself. Candidates the network calls noise are dropped **before** the transform, so only accepted cells reach atlas space and the regional counts.

The result is cached next to the registration folder as `<name>_classification.mat`:

| Field | Contents |
| :--- | :--- |
| `isgood` | `[Ncells × 1]` logical — the cells that were kept |
| `labels`, `scores` | raw network output per candidate |
| `classnames`, `goodclass` | the network's classes and which one means "cell" |
| `fracgood`, `ncells` | fraction kept, and how many were classified |
| `netname`, `netpath`, `netaccuracy` | which network produced this |

**A later run reuses that file instead of re-running the network.** Pass `'reclassify', true` to force a re-run after retraining; a cache that no longer matches the detection count is discarded automatically. The `_atlas.mat` output records the classification it was filtered with.

### Notes and limits

* **Classification needs cell images**, so it applies to `.mat` detections saved with `savecellimages = true`. CSV and XML point sets carry no images and are transformed unfiltered, with a warning.
* **Candidates are classified in chunks** (10,000 at a time) so a million-cell brain fits in memory. The image scaling is measured once over the whole file and reused per chunk, so the result is identical to classifying in one pass.
* **Retrain when the staining changes.** The network learns *your* artifacts and will not transfer perfectly to another label or microscope. Relabelling a few hundred candidates takes under an hour.
* **Check what it rejects, not just what it keeps** — `visualizeNetworkClassification` on your validation set is the fastest way to spot a network that learned the wrong cue.

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
│                                            #   (+ cell_images if savecellimages=true)
├── chan_1_cell_locations_labeled.mat        # Your CellLabelingTool labels
├── chan_1_cell_locations_classification.mat # Cached CNN verdict (if 'network' given)
├── chan_1_cell_detections/                  # Debug images (if debug=true)
├── dim1/2/3_initial_registration.png        # Initial alignment check
├── *_affine_registration.png                # After affine step
├── *_bspline_registration.png               # After B-spline step
└── volume_registered/
    ├── chan01_intensities.mat               # Regional intensity & volume stats
    ├── chan01_intensities.csv               # (if writetocsv=true)
    ├── chan01_cellcounts.mat                # Regional cell counts
    ├── chan_1_cell_locations_atlas.mat      # Cells in atlas coordinates
    └── cell_counts_by_region.csv           # Regional cell counts (if writetocsv=true)
```

The trained networks themselves live wherever you point `netowrksavepath` in `demos/trainClassificationNetwork.m`; `<date>_CellClassifierNet.mat` files are excluded from version control.
