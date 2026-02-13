# Installation

LightSuite relies on external tools and repositories. To install LightSuite and stay up-to-date, close the LightSuite repository and add it to MATLAB's path. If you're not familiar with GitHub, the easiest starting point is using [GitHub Desktop](https://github.com/apps/desktop). The cleanest way to add a repository to MATLAB's path is by adding the following line to `~\Documents\MATLAB\startup.m`:

```addpath(genpath('~\Documents\GitHub\LightSuite'))```

Follow these steps to set up your environment correctly: 

### 1. MATLAB Requirements
You need **MATLAB R2022b** or newer. Ensure the following toolboxes are installed:

* Computer Vision Toolbox
* Image Processing Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox
* [Bioformats Image Toolbox](https://ch.mathworks.com/matlabcentral/fileexchange/129249-bioformats-image-toolbox)
* [Bioformats functions](https://downloads.openmicroscopy.org/bio-formats/5.5.0/artifacts/bfmatlab.zip)

### 2. Install repository dependencies
Clone the following two repositories and add them to MATLAB's path: 

* [matlab_elastix](https://github.com/dimokaramanlis/matlab_elastix), updated fork based on Rob Campbell's repository
* [bcpd](https://github.com/ohirose/bcpd), by Osamu Hirose

### 3. Install Elastix
LightSuite calls Elastix from the command line, so its executables must be globally accessible. 

1. Download the [version 5.1.0 executables](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) and unzip the archive to a permanent location on your computer (e.g., C:\elastix or /opt/elastix).

2. Add Elastix to the system PATH: This is the most critical step. You must add the directory containing the elastix and transformix executables to your system's PATH environment variable.
	- Windows: Search for "Edit the system environment variables," click "Environment Variables," and add a new System variable named elastix with the path to your bin directory (e.g., C:\elastix).
	- Linux/macOS: Add export PATH="/path/to/your/elastix/bin:$PATH" to your shell configuration file.
	
3. Verification: To confirm that Elastix is installed correctly, open a new terminal or command prompt and type elastix --version. You should see the version information printed to the screen. If you get a "command not found" error, the PATH is not set correctly.


### 4. Download Atlas Data

LightSuite uses the [2020 version](https://alleninstitute.github.io/abc_atlas_access/descriptions/Allen-CCF-2020.html) of the Allen Brain Atlas for brain registration. To be able to use the software you have to download all [volume files](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#image_volumes/Allen-CCF-2020/20230630/) along with all relevant [metadata](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/Allen-CCF-2020/20230630/) and make them available to the MATLAB path.

For spinal cord registration, LightSuite uses the [Fiederling et al 2021](https://www.sciencedirect.com/science/article/pii/S2667237521001260) spinal cord atlas. Download the compressed file from [here](https://data.mendeley.com/datasets/4rrggzv5d5/1), extract its contents and make them available to the MATLAB path.

### 5. Check installation
Finally, run `check_lightsuite_installation.m` to verify whether you have installed everything correctly. If you encounter no errors, you're good to go.