LightSuite is a comprehensive, semi-automated software pipeline designed for the end-to-end analysis of whole-brain imaging data acquired via light-sheet and widefield slice microscopy. The pipeline addresses the significant challenge of processing and quantifying cellular-level data from datasets in the 100GB scales. Its core functionalities are twofold: first, the robust registration of raw experimental image volumes to the standardized Allen Common Coordinate Framework (CCF) v3 for the mouse brain; and second, the automated 3D detection and quantification of labeled cells within this standardized anatomical context. LightSuite enables researchers to transition from raw image stacks to quantitative, region-specific cell counts, providing a crucial tool for systems neuroscience.

## Registration

![Example bspline registration](./images/example_bspline.png)

## MATLAB requirements

Cell counting is done on the GPU, needs to improve.

Requires:
	- MATLAB >= R2022b
	- Computer Vision Toolbox
	- Image Processing Toolbox
	- Optimization Toolbox
	- Parallel Computing Toolbox
	- Statistics and Machine Learning Toolbox

## elastix installation

LightSuite calls Elastix from the command line, so its executables must be globally accessible. 

1. Download the [version 5.1.0 executables](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) and unzip the archive to a permanent location on your computer (e.g., C:\elastix or /opt/elastix).

2. Add Elastix to the system PATH: This is the most critical step. You must add the directory containing the elastix and transformix executables to your system's PATH environment variable.
	- Windows: Search for "Edit the system environment variables," click "Environment Variables," and add a new System variable named elastix with the path to your bin directory (e.g., C:\elastix).
	- Linux/macOS: Add export PATH="/path/to/your/elastix/bin:$PATH" to your shell configuration file.
	
3. Verification: To confirm that Elastix is installed correctly, open a new terminal or command prompt and type elastix --version. You should see the version information printed to the screen. If you get a "command not found" error, the PATH is not set correctly.



## Allen Atlas
[2020 version](https://alleninstitute.github.io/abc_atlas_access/descriptions/Allen-CCF-2020.html), download all files [here](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#image_volumes/Allen-CCF-2020/20230630/)
and all metadate [here](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/Allen-CCF-2020/20230630/)


