# Spinal Cord Volume Registration

This module processes lightsheet spinal cord datasets and registers them to the Fiederling et al. spinal cord atlas. The main challenge specific to spinal cord data is that the tissue is curved — this module handles straightening the cord before atlas registration.

## Before You Start

You will need:
* A stitched lightsheet TIFF dataset of the spinal cord.
* [Elastix](https://elastix.lumc.nl/) installed and on your system path.
* The spinal cord atlas files (loaded via `loadSpinalCordAtlas()`).

---

## Workflow Overview

The analysis is driven by the main script: `ls_analyze_spinal_cord`.

| Step | Stage | Type |
| :---: | :--- | :--- |
| 1 | Load and prepare sample | Automated |
| 2 | Manual spinal cord straightening (Spinal Cord Aligner GUI) | **Manual** |
| 3 | Initial registration | Automated |
| 4 | Manual alignment refinement (Control Point GUI) | **Manual** |
| 5 | Deformable registration | Automated |
| 6 | Volume registration to atlas | Automated |

---

## 1. Load and Prepare Sample

**Functions:** `readSpinalCordSample(dpspinesample, sampleres)` → `prepareCordSampleForRegistration(cordvol, opts)`

These functions load your raw TIFF data at the specified resolution and prepare a downsampled registration volume, similar to the lightsheet brain pipeline. The registration options are saved to `regopts.mat`.

---

## 2. Manual Spinal Cord Straightening

**Function:** `spinal_cord_aligner(regopts)`

Before registration, the curved spinal cord must be "straightened" to match the linear coordinate system of the atlas. This GUI lets you trace anatomical landmarks along the cord, and LightSuite fits a smooth spline through your markings to model the cord's shape.

### The Interface

The GUI has three panels:

1. **Left (Image):** The transverse (cross-sectional) view of the current spinal cord slice.
2. **Middle (Position Plot):** Tracks the X/Y center position of your markers across all slices.
3. **Right (Angle Plot):** Tracks the rotation angle (θ) of the cord at each slice.

Use the side plots to verify that your markings form a smooth curve. A sudden spike in either plot usually means a misplaced point.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **Left Click** | Mark **Anterior/Front** (green point) |
| **Right Click** | Mark **Posterior/Back** (red point) |
| **Middle Click** | Mark **Center Canal** (blue point) |
| **← / →** or **Scroll Wheel** | Navigate slices |
| **Spacebar** | Jump to the largest gap in your annotations |
| **p** | Toggle prediction visibility (cyan crosses/lines) |
| **c** | Clear all points on the current slice |
| **x** | Delete the single point nearest to the cursor |
| **l** | Set regularization parameters (spline stiffness) |
| **s** | Save alignment data |
| **Click on side plot** | Jump directly to that slice |

### Workflow

1. Mark a few slices near the start of the cord (anterior end).
2. Press **Spacebar** to jump to the largest unmarked gap and mark that slice.
3. Repeat until the cyan prediction lines track the cord smoothly throughout the volume.
4. Check the side plots — the position and angle curves should be smooth and continuous.
5. Press **s** to save.

**Tips:**
* You do not need to mark every slice. A sparse set of evenly distributed annotations is sufficient — the spline fills in the gaps.
* If you see a sharp spike in the position or angle plot, you have a misplaced point. Use **x** to delete it and remark that slice.
* Use **l** to adjust the regularization (`lambda_pos`, `lambda_ang`) if the spline is too stiff or too flexible. Larger values produce a smoother fit; smaller values follow your points more closely.

**Output:** `spinal_alignment_opt.mat` — contains the fitted spline parameters (`fit_x`, `fit_y`, `fit_theta`) and your raw user clicks.

---

## 3. Initial Registration

**Function:** `initializeCordRegistration(regopts)`

After straightening, LightSuite automatically computes an initial coarse alignment between the straightened cord volume and the spinal cord atlas using a point-cloud-based similarity transform. This gives the manual control point step a sensible starting point.

---

## 4. Manual Alignment Refinement (Control Point GUI)

**Function:** `matchControlPointsSpine(regopts)` — internally calls `matchControlPoints_unified(opts)`

This GUI lets you refine the alignment by clicking matching anatomical landmarks in the sample and atlas side by side.

### Interface

* **Left panel (80%):** Your histology/sample slice.
* **Right panel (80%):** The corresponding atlas slice.
* **Side panels (20%):** Maximum intensity projections (MIPs) of the sample and atlas along the long axis, with a line cursor showing the current slice position.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **Click (Left panel)** | Place a control point on the sample |
| **Click (Right panel)** | Place the corresponding point on the atlas |
| **← / →** | Navigate to the previous / next slice |
| **Scroll Wheel** | Scroll atlas slice plane independently |
| **Backspace** | Delete the last point added |
| **c** | Clear all points on the current slice |
| **Spacebar** | Toggle the red atlas boundary overlay |
| **Enter** | Jump to a specific slice number |
| **s** | Save and exit |

### Tips for Good Results

* **Minimum 16 matched pairs** are required for a robust 3D affine fit. Place points across the full length of the cord (rostral to caudal) and in different radial positions (dorsal, ventral, lateral).
* The fit quality (MSE) updates live in the top banner — a lower MSE indicates better correspondence.
* **Good landmarks** include the central canal, the dorsal horn boundaries, and white matter funiculi outlines.
* Use the **Scroll Wheel** on the atlas panel if the correct anatomical plane isn't shown — the atlas and sample can be scrolled independently.
* If a tissue region is damaged or missing, estimate where it *should* be and place a point there anyway. This prevents the deformable registration from bending the atlas into the damaged void.

**Output:** `corresponding_points.mat` — 3D affine transform and all control point arrays.

---

## 5. Deformable Registration

**Function:** `multiobjCordRegistration(regopts, control_point_weight)`

Refines the alignment using Elastix B-spline registration, incorporating both your manual control points and the initial automatic landmarks. Both forward (atlas → sample) and reverse (sample → atlas) transforms are computed.

**Outputs:**
* `transform_params.mat` — complete registration parameters
* Elastix transform files (forward and reverse)
* Diagnostic registration images

---

## 6. Volume Registration to Atlas

**Function:** `generateRegisteredCordVolume(regopts, transform_params)`

Warps all channels into atlas space and computes regional statistics per anatomical area.

### Atlas Regions

The spinal cord atlas is organized into two levels:

* **Division level:** Gray Matter (GM) and White Matter (WM).
* **Structure level:** Laminas I–X and the dorsal (`df`), lateral (`lf`), and ventral (`vf`) funiculi.

Statistics are aggregated along the rostrocaudal axis in anatomically defined segments.

**Output files** (in `volume_registered/`):
* `chan0X_intensities.mat` — median signal intensity and volume per region and segment
* Registered volume slices

---

## Post-Registration Analysis

After registration, use `reorganizeSpinalCordAreas()` to aggregate the regional data to your desired anatomical level:

```matlab
% Aggregate to Gray Matter / White Matter level
resout = reorganizeSpinalCordAreas(signals, volumes, parcelinfo, avinds, 'division');

% Aggregate to Lamina / Funiculus level
resout = reorganizeSpinalCordAreas(signals, volumes, parcelinfo, avinds, 'structure');
```

The `resout` struct contains `counts`, `volumes`, `signal`, `names`, and `indices` fields for downstream plotting and statistics. See `demos/example_analysis_spinal_cord.m` for a complete worked example including heatmaps and intensity profile plots.

---

## Output File Summary

```
<lsfolder>/
├── regopts.mat                        # Registration options (updated at each stage)
├── spinal_alignment_opt.mat           # Straightening spline parameters
├── corresponding_points.mat           # Manual control points and initial affine
├── transform_params.mat               # Final registration parameters
├── elastix_forward/                   # Forward B-spline transforms
├── elastix_reverse/                   # Reverse B-spline transforms
└── volume_registered/
    ├── chan0X_intensities.mat         # Regional intensity and volume stats
    └── [registered volume slices]
```
