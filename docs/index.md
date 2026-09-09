# Welcome to LightSuite

**LightSuite** is an open-source MATLAB framework for **scalable atlas registration and cell detection across the central nervous system**. It combines damage-resilient 3D point-cloud registration, user-steerable non-rigid registration, and GPU-accelerated local background normalization for fast, illumination-invariant cell detection — on a standard workstation, at cohort scale.

## What can I do with LightSuite?

| Modality | What the module does |
| :--- | :--- |
| [**Whole-brain volumes**](usage_lightsheet_brain.md) | Detect cells at native resolution and register the volume to the Allen CCFv3; export per-region counts, densities, soma diameters and intensities per hemisphere. Any volumetric modality — light-sheet, fMOST, serial two-photon. |
| [**Spinal cord**](usage_spinal_cord.md) | Computationally straighten and untwist a curved cleared cord, then register it to the Fiederling et al. (2021) atlas. |
| [**Wide-field slices**](usage_slice.md) | Order, align and non-rigidly register serial sections, then assemble them into a 3D volume in atlas space with 2D cell detection. |
| [**Functional ultrasound (fUS)**](usage_fusi.md) | Align repeated sessions to a within-mouse seed, register that to a vascular CCFv3 template, and carry activation maps and timecourses into atlas space — or onto the cortical flatmap. |
| [**Probe & implant tracing**](usage_tracing.md) | Localize optical fibers / GRIN lenses and Neuropixels probes in atlas space, with the regions and fluorescence beneath the implant. |

Registration does not need a dedicated autofluorescence channel: a single reporter-fluorescence volume usually carries enough anatomical contrast, halving acquisition time and storage.

## Hardware Requirements

A **dedicated GPU** is recommended for large 3D volumes (band-pass filtering and cell detection). The **slice module** runs comfortably without one. Reference timings on an i9-10900X / 64 GB / RTX A4000: whole-brain cell detection 3–4 h, automated registration 5–6 min, landmark curation 10–20 min.

## Supported Data Formats

| Modality | Accepted input |
| :--- | :--- |
| Whole-brain volumes | Single-channel 2D TIFF planes, multi-channel volume TIFFs, or per-channel volume TIFFs split across files. Acquisition modality does not matter — light-sheet, fMOST, serial two-photon and comparable stacks all work |
| Spinal cord | Low-resolution whole-cord volumes; channels together in one TIFF or separate |
| Slices | 2D TIFF planes (one per slice), or AxioScan **`.czi`** output |
| fUS | Repeated 3D power-Doppler scans as MATLAB arrays (a `*_FUS.mat` reader is provided; any loader works) |

Support for other brain orientations is planned.

## Getting Started

1. **Install** MATLAB dependencies and external tools (Elastix) — see [Installation](installation.md).
2. **Skim** [How it works](how_it_works.md) for the registration and detection concepts every workflow builds on.
3. **Configure** by editing the `opts` struct at the top of the relevant demo script (cell diameter, paths, channels).
4. **Pick your workflow** from the table above.
