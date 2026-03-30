# Slice Analysis

The Slice Module is optimized for conventional wide-field microscope data where you have individual coronal slices rather than a continuous 3D volume. 

The primary pipeline is executed via the `ls_analyze_slice_volume` script. Unlike lightsheet data, slice registration inherently requires manual curation due to missing tissue, arbitrary cutting angles, and uneven spacing.

## Pipeline Workflow & Algorithmic Steps

The pipeline alternates between automated algorithmic processing and manual user interventions. 

### 1. Volume Initialization & Curation
* **Automated Stacking:** The pipeline first reads individual 2D images (e.g., `.czi` files) and stacks them into a preliminary 3D slice volume.
* **Manual Curation (`SliceOrderEditor`):** A GUI is presented to allow the user to reorder slices, flip them horizontally (if mounted upside down), and discard damaged or irrelevant slices. A new, curated volume is then generated.

### 2. Initial Coarse 3D Alignment (`alignSliceVolume`)
Before 2D slice-to-slice registration, the entire ordered stack is coarsely aligned to the Allen Brain Atlas in 3D.
* **Signal Extraction:** The volume is median-filtered to estimate background, and a $\Delta F/F$ (signal-over-background) image is computed. 
* **Point Cloud Generation:** High-intensity pixels are thresholded (top 1%) and extracted as a 3D point cloud representing the sample. A corresponding point cloud is generated from the Allen Atlas template.
* **Iterative Point Cloud Registration:** The sample point cloud is registered to the atlas point cloud using sequential optimization steps, beginning with rigid transformations and followed by affine transformations.

### 3. Cutting Angle Correction & Control Points (Manual)
* **Angle Estimation:** The user uses a GUI (`determineCuttingAngleGUI`) to correct the slicing angle relative to the perfect coronal plane of the atlas.
* **Control Point Matching:** The user maps specific anatomical landmarks in the 2D sample slices to the corresponding 2D slices in the Allen Atlas. This anchors the subsequent elastic registration.

### 4. Refined 2D Slice-to-Atlas Registration (`registerSlicesToAtlas`)
Each slice is individually registered to its corresponding AP (Anterior-Posterior) plane in the atlas.
* **Interpolation:** If control points are missing for certain slices, the algorithm interpolates the expected atlas indices using the available user-provided points.
* **Affine Registration:** For each slice, an initial 2D affine transformation is calculated. If sufficient control points (>4) are provided, a geometric transform is fitted. If absent, an image-based point-cloud affine fit is used.
* **Elastix B-Spline Registration:** The affine-transformed atlas image is further deformed to match the histological slice using a non-linear B-spline registration powered by Elastix (`bsplineRegisterSlice`), heavily weighting the user's control points. The transformation is then inverted to map sample coordinates back to atlas space.

### 5. 2D Cell Detection (`extractCellsFromSliceVolume`)
Cell detection is performed on a per-slice basis.
* **Background Subtraction:** A large median filter (based on expected cell diameter) calculates local background, which is used to generate a normalized $\Delta F/F$ intensity image.
* **Blob Detection:** A 2D detection algorithm (`cellDetector2D`) identifies cell centroids based on user-defined expected diameter and Signal-to-Background (SBR) thresholds.

### 6. Mapping to Atlas Space
Finally, the detected 2D cell coordinates are mapped backwards through the previously calculated Elastix B-spline and Affine transforms into the standard 3D Allen Atlas coordinate space (`slicePointsToAtlas`), allowing for regional cell counting.