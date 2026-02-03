# Installation

LightSuite is a MATLAB-based pipeline that relies on external tools and repositories. Follow these steps to set up your environment correctly.

## 1. MATLAB Requirements
You need **MATLAB R2022b** or newer. Ensure the following toolboxes are installed:

* Computer Vision Toolbox
* Image Processing Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

## 2. External Dependencies
LightSuite requires three external repositories to function. You should clone these into a folder where MATLAB can access them (e.g., your `Documents/MATLAB` folder).

### A. MATLAB Helper Packages
Clone the following two packages (originally by Rob Campbell) which manage Elastix communication and YAML parsing:

1.  **matlab_elastix**:
    ```bash
    git clone [https://github.com/raacampbell/matlab_elastix.git](https://github.com/raacampbell/matlab_elastix.git)
    ```
2.  **yamlmatlab**:
    ```bash
    git clone [https://github.com/raacampbell/yamlmatlab.git](https://github.com/raacampbell/yamlmatlab.git)
    ```

### B. BCPD (For Spinal Cord)
If you plan to use the Spinal Cord module, you must install **Bayesian Coherent Point Drift (BCPD)**.

1.  **Download/Clone**:
    ```bash
    git clone [https://github.com/ohirose/bcpd.git](https://github.com/ohirose/bcpd.git)
    ```
2.  **Binaries**:
    * **Windows**: The repository usually contains a `win/bcpd.exe` folder. Ensure this executable is accessible.
    * **Linux/Mac**: You may need to compile the source code using `make` following the instructions in the BCPD repository.

!!! important "Update MATLAB Path"
    Once downloaded, add `matlab_elastix`, `yamlmatlab`, and the `bcpd` folder to your **MATLAB Path** (Home -> Set Path -> Add with Subfolders).

## 3. Install Elastix
LightSuite uses **Elastix** for image registration.

1.  **Download**: Get the **5.1.0** version executables from [GitHub](https://github.com/SuperElastix/elastix/releases).
2.  **Unzip**: Extract to a permanent location (e.g., `C:\elastix`).
3.  **Add to PATH (Critical)**:
    * **Windows**: Add the `elastix\bin` folder to your System Environment Variable `Path`.
    * **Linux/macOS**: Add `export PATH="/path/to/elastix/bin:$PATH"` to your shell config.

## 4. Atlas Data
LightSuite requires reference atlas files to perform registration.

### Allen Brain Atlas (CCFv3)
For whole-brain lightsheet or slice analysis:
1.  Download the **Allen CCF v3 (2020)** volume files and metadata.
2.  Place them in a dedicated folder (e.g., `D:\AllenAtlas\`).

### Spinal Cord Atlas
For the spinal cord module, you need the specific spinal cord reference volume:
1.  Obtain the **Allen Cord Atlas** (e.g., `allen_cord_20um_v1.1`).
2.  Place this in your atlas directory (e.g., `D:\AllenAtlas\extra_spine\`).

## 5. Install LightSuite
Finally, clone the main repository.