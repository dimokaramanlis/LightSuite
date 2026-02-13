# Spinal cord volume registration


## Manual spinal cord straightening

This module handles the specific challenges of spinal cord data, including straightening curved samples and registering them to the spinal cord atlas.

## 1. Manual Spinal Cord Straightening

Before registration, the curved spinal cord volume must be "straightened" to match the linear coordinate system of the atlas. This is achieved using the **Spinal Cord Aligner GUI**.

### Launching the Tool
Run the straightening function (typically called via `spinal_cord_aligner(opts)`). You must provide a struct containing the registration volume and save folder.

### The Interface
The GUI consists of three main panels:
1.  **Left (Image):** The transverse view of the spinal cord slice.
2.  **Middle (Position Plot):** Tracks the X/Y position of your markers across slices.
3.  **Right (Angle Plot):** Tracks the rotation angle ($\theta$) of the cord.

### Workflow
Your goal is to mark anatomical landmarks on various slices. The algorithm fits a smooth spline through your points to predict the shape of the entire cord.

**Controls:**

| Key / Action | Function |
| :--- | :--- |
| **Left Click** | Mark **Anterior/Front** (Green point) |
| **Right Click** | Mark **Posterior/Back** (Red point) |
| **Middle Click** | Mark **Center Canal** (Blue point) |
| **Space** | Jump to the largest gap in your annotations (auto-navigation) |
| **L** | Set Regularization parameters (stiffness of the fit) |
| **P** | Toggle prediction visibility (Cyan crosses/lines) |
| **C** | Clear all points on the *current* slice |
| **X** | Delete the single point nearest to the cursor |
| **S** | Save the alignment data |

**Tips:**
* You do not need to mark every slice. Mark a slice, press `Space` to jump to a gap, and mark again. The spline interpolation will fill in the rest.
* Use the side plots to verify smoothness. If you see a sharp spike in the Position or Angle plots, you likely have a misplaced point.

## Manual Alignment Correction (Control Point GUI)

The **Manual Alignment GUI** is a tool within LightSuite designed to refine the spatial registration between your experimental histology data and the reference atlas (CCF). While the automated pipeline handles most samples, brains with significant tissue damage, clearing artifacts, or non-linear deformations often require manual intervention to achieve a high-fidelity fit.

### Overview

The interface displays your experimental data alongside the reference atlas. By manually identifying corresponding anatomical landmarks in both views, you update the affine transformation matrix. This "globally" aligns the brain before any subsequent non-linear warping steps.

* **Left Panel (Sample):** Displays your histology/sample slice.
* **Right Panel (Atlas):** Displays the corresponding slice from the reference atlas (CCF).
* **Top Banner:** Displays controls and the **Current Fit Status** (MSE score and number of control points).

### Workflow

1.  **Launch the GUI:** Run the `matchControlPointsSpine` function with your registration options.
2.  **Navigate Slices:** Use `Left Arrow` / `Right Arrow` to scroll through your sample slices.
    * **Note on Randomization:** You will notice that the slices appear in a **randomized order**, rather than sequentially (e.g., jumping from anterior to posterior). This is intentional. It encourages you to place control points evenly throughout the entire volume, preventing you from clustering all points in a single region or "over-fitting" one specific area.
3.  **Find a Landmark:** Identify a distinct anatomical feature visible in both the Sample and the Atlas (e.g., a specific ventricle shape, fiber tract, or nuclei boundary).
4.  **Place Control Points:**
    * **Click** the feature on the **Sample** (Left) image.
    * **Click** the corresponding location on the **Atlas** (Right) image.
    * *Points are color-coded and numbered. A matched pair shares the same number.*
5.  **Refine & Repeat:**
    * Add points across different slices to constrain the 3D volume.
    * **Minimum Points:** The fit calculates automatically once **16 points** are placed.
    * **Check Quality:** Watch the **MSE (Mean Squared Error)** in the top banner. A lower MSE indicates a tighter correspondence between your points.
6.  **Verify Alignment:** Press `Space` to toggle the red Atlas overlay on/off to visually verify the improved fit.
7.  **Save:** Press `S` to save the new alignment transformation.

### Controls Reference

| Key / Action | Function |
| :--- | :--- |
| **Left Click** | Add a control point (Sample or Atlas). |
| **← / →** | Previous / Next Sample slice (Randomized order). |
| **Scroll Wheel** | Scroll through **Atlas** slices (independent of sample). |
| **Backspace** | Delete the last added point. |
| **C** | Clear all points on the *current* slice. |
| **Space** | Toggle the red Atlas boundary overlay. |
| **Enter** | Jump to a specific slice index. |
| **S** | **Save** the alignment and export points. |

### Tips for Best Results

* **Handling Damaged Tissue:** If a region is damaged or missing in your sample, **do not ignore it**.
    * Try to "imagine" where the anatomical structure *should* be if the tissue were intact.
    * Place a point at this estimated location in the sample, and the corresponding true location in the atlas.
    * *Why?* This anchors the transformation and prevents the subsequent non-linear registration from "bending" the atlas into the damaged void, ensuring the atlas maintains a realistic shape.
* **Spread your points:** Do not cluster points in one area. The randomized slice navigation helps with this, but you should still consciously aim to mark landmarks in the center, periphery, dorsal, and ventral regions.
* **Atlas Scrolling:** If the correct anatomical plane isn't visible on the right, use the **Scroll Wheel** to move the Atlas slice forward/backward until it matches the Sample's Z-position.