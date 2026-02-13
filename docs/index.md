# Welcome to LightSuite

**LightSuite** is a comprehensive MATLAB-based pipeline designed for the registration and analysis of large-scale microscopy datasets. It provides modular workflows for whole-brain lightsheet volumes, coronal slice series, and spinal cord data, bridging the gap between raw microscopy images and anatomical reference atlases (CCF).

## What can I do with LightSuite?

LightSuite automates the complex tasks of mapping experimental data to standard anatomical coordinates and quantifying labeled cells.

* **Whole-Brain Lightsheet Analysis:** Process continuous 3D volumes. The pipeline handles preprocessing (median filtering, binary conversion), automated cell detection (SNR-based local maxima), and registration to the Allen Brain Atlas.
* **Spinal Cord Analysis:** Specialized tools for straightening and registering spinal cord volumes. It includes a dedicated GUI for defining the central canal and anterior/posterior axes to unroll and map the cord before registration.
* **Slice Analysis:** Optimized for conventional wide-field microscope data (e.g., coronal slices). It registers individual 2D planes to the atlas and outputs registered image stacks and cell coordinates.

## Hardware Requirements

* **Standard Workstations:** The **Slice Analysis** module is optimized for efficiency and has been tested on standard computers without GPUs.
* **High-Performance Workstations:** For **Large-scale lightsheet volumes**, we recommend a system with a dedicated GPU to accelerate 3D operations (such as spatial band-pass filtering and cell detection).

## Supported Data Formats

### 1. Large-scale Lightsheet Volumes
The pipeline currently supports data formatted as a series of **2D TIFF planes** representing a single channel, sliced axially. Support for other brain orientations is planned.

### 2. Spinal Cord Data
We support low-resolution whole cord volumes. Channels can be stored within the same TIFF volume or separated.

### 3. Slice Volumes
We support:
* A series of **2D TIFF planes** (one file per slice).
* Direct output from AxioScan scanners (**`.czi`** files).

## Getting Started

1.  **Installation:** Follow the instructions in [Installation](installation.md) to set up MATLAB dependencies and external tools (Elastix).
2.  **Configuration:** LightSuite uses script-based configuration. You will adjust parameters (such as cell diameter or file paths) directly within the analysis scripts.
3.  **Select your Workflow:**
    * [Lightsheet Brain Analysis](usage_lightsheet_brain.md)
    * [Spinal Cord Analysis](usage_spinal_cord.md)
    * [Slice Analysis](usage_slice.md)