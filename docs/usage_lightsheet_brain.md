# Lightsheet Whole-Brain Analysis

This module is designed for processing whole-brain datasets acquired via lightsheet microscopy. It handles large-scale volumetric data, performing artifact removal, cell detection, and registration to the Allen Brain Atlas (CCF).

## Workflow Overview

The analysis is primarily driven by the main script: `ls_analyze_lightsheet_volume`.

The pipeline consists of three distinct stages:
1.  **Preprocessing:** Data conversion and artifact removal.
2.  **Cell Detection:** SNR-based identification of cells in 3D.
3.  **Registration:** Aligning the sample to the atlas.

---

## 1. Preprocessing

The preprocessing stage prepares raw microscope data for efficient analysis. This is handled by the `preprocessColmVolumeNew` function.

### Key Operations
* **Artifact Removal:** A 3x3 median filter is applied to each slice to remove salt-and-pepper noise and imaging artifacts.
* **Data Conversion:** The filtered volume is written to a large binary file. This allows for faster random access during the cell detection stage compared to reading individual TIFF files.
* **Downsampling for Registration:** A lower-resolution version of the volume (e.g., `sample_register_Xum.tif`) is generated and saved specifically for the registration step. This volume is intensity-rescaled to match the bit-depth requirements of the registration tools.

**Configuration:**
Parameters such as pixel size (`opts.pxsize`) and registration resolution (`opts.registres`) are defined directly in the main analysis script options.

---

## 2. Cell Detection

Cell detection is performed in batches on the high-resolution data using the `cellDetectorNew` function.

### Algorithm
The detection pipeline uses a 3D signal-to-noise ratio (SNR) approach:
1.  **3D Bandpass Filtering:** The volume is filtered to enhance objects of a specific size while suppressing high-frequency noise and low-frequency background.
2.  **Robust Noise Estimation:** The background standard deviation is estimated using a robust statistical method to define the noise floor.
3.  **Local Maxima Detection:**  Candidate cells are identified as local maxima that exceed a specific SNR threshold.
4.  **Morphological Filtering:** Candidates are further filtered based on shape (ellipticity) and size constraints to exclude artifacts.

### Key Parameters
* **`avgcellradius`**: The expected radius of the cells in pixels. This is the most critical parameter for the bandpass filter.
* **`sdthres`**: The SNR threshold. Lower values detect dimmer cells but may increase false positives.

---

## 3. Registration

Registration maps your experimental volume to the Allen Brain Atlas (CCF).

* **Input:** Uses the downsampled volume generated during preprocessing.
* **Method:** The pipeline supports both fully automated registration and a manual intervention mode.
* **Manual Adjustment:** While automation is possible, manual adjustment is highly recommended to ensure high-fidelity mapping. The user can visually inspect the alignment and adjust parameters if the automated fit is insufficient.