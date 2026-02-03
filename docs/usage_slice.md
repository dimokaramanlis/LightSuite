# Slice Analysis

The Slice Module is optimized for conventional wide-field microscope data where you have individual coronal rather than a continuous 3D volume.

## Workflow

1.  Run the slice analysis script:
    ```matlab
    ls_analyze_slice_volume
    ```
2.  **Manual Intervention**: Unlike the fully automated lightsheet pipeline, slice registration often requires manual matching.
    * The GUI will present each slice.
    * You will need to confirm the approximate Anterior-Posterior (AP) position of the slice in the Atlas.
3.  **Registration**: The script calculates the warp from the 2D slice to the corresponding 2D atlas plane.
