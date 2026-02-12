# Manual Alignment Correction (Control Point GUI)

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