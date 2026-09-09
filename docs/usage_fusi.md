# Functional Ultrasound (fUS)

LightSuite registers **power-Doppler functional ultrasound (fUS)** volumes into the Allen CCFv3 and carries functional results — activation maps, per-area timecourses, whole recordings — into atlas space. Three things make fUS harder than fixed tissue, and the module addresses each:

* **The contrast is vascular, not cytoarchitectonic.** Registering power-Doppler against the Allen average template compares two unrelated images. LightSuite instead registers against a **vascular template** (Brunner et al., 2021) that has itself been registered into CCFv3 and inherited the Allen annotation, so image similarity is computed between vessels and vessels.
* **The probe is positioned freehand.** Each session images a different slab, so a mouse has no single anatomical volume until you build one. Sessions are rigidly aligned to a within-mouse **seed** scan and averaged into an anatomy scan defining a common per-mouse space.
* **Resolution is coarse** (typically ≥100 µm), so the non-rigid fit is constrained to a **1.6 mm isotropic** B-spline grid to prevent overfitting to vascular fluctuations, and guided by user landmarks — a median of 219 per mouse in the paper.

Everything else is shared with the rest of LightSuite: the same [control-point GUI](#3-place-matched-control-points), the same [similarity → affine → B-spline transform chain](how_it_works.md#registration-from-sample-to-atlas), the same `regopts.mat` / `transform_params.mat` on disk.

Against an affine baseline optimized within LightSuite on the same landmarks, the non-rigid fit suppresses activation leakage across atlas borders, raises bilateral response correlations, and raises pairwise inter-animal correlations for both vascular anatomy and activation maps.

**Demo script:** `demos/ls_analyze_fusi.m` — a runnable, annotated version of everything on this page.

---

## Before You Start

You will need:

* **Elastix 5.1.0** on the system `PATH` and the MATLAB dependencies from the [installation guide](installation.md).
* **The Allen CCF (2020)** on the MATLAB path (`annotation_10.nii.gz`, `parcellation_to_parcellation_term_membership.csv`).
* **A vascular fUS template** — loaded by `loadAtlasInfo('allen2020fusi_50um')`. `src/fusi/prepare_fusi_atlas.m` builds one: the high-resolution power-Doppler template (50 µm isotropic) is registered into CCFv3 with the same point-cloud-initialized, multi-metric B-spline pipeline used for cleared brains, the Allen annotation is transferred onto the warped template, and the result is symmetrized across the midline to remove residual left-right asymmetry. Point `loadAtlasInfo` at wherever you keep it.
* **(Optional, flatmaps only)** a Python environment with `allensdk` and `ccf_streamlines`, plus the Allen flatmap resources (`flatmap_butterfly.h5` / `.nrrd`, `surface_paths_10_v3.h5`, `labelDescription_ITKSNAPColor.txt`, `manifest.json`).

### Data conventions

| | |
| :--- | :--- |
| **Raw session volume** | `[d1 d2 d3]`, coronal slices along **dimension 3** |
| **Raw voxel size** | anisotropic, in **mm**, e.g. `[0.1971 0.15 0.15]` |
| **Functional recording** | `[d1 d2 d3 nframes]` in the same geometry, time last |
| **Anatomy / seed space** | isotropic at `opts.atlas_res` (50 µm), still in the **native session orientation** |
| **Atlas space** | CCFv3, `AP × DV × ML`, 50 µm isotropic — `[264 160 228]` |

The sample→atlas axis permutation is chosen once per mouse and stored; nothing before that step is reoriented.

---

## Pipeline Overview

```
per-session volumes ──(1) buildFusiAnatomy──▶ anatomy / seed space
                                                    │
                                    (2) setupFusiOptions  (orientation)
                                    (3) matchControlPoints_minimal  (landmarks)
                                    (4) multiobjRegistrationFusi    (affine + B-spline)
                                                    │
                                                    ▼
functional recording ──(5a) fusiRecordingToAnatomy──▶ seed space
                       ──(5b) applyFusiTransforms ───▶ ATLAS SPACE
                                                    │
                                    (6) fusiVolumeToFlatmap ──▶ cortical flatmap
```

Steps 1–4 run **once per mouse** and cache their results; re-running the script walks past the GUIs without asking again. Step 5 runs once per recording.

---

## 1. Rigid Alignment Across Sessions

**Function:** `buildFusiAnatomy(opts, sessionvols, voxelsize_mm, sessionnames)`

Each session first has to be reduced to one 3-D volume — a robust temporal average, or *session template*. `loadFusiSessionTemplate` does this for recordings stored as a `*_FUS.mat` with an `I = [nvox × nframes]` array: ten candidate volumes are formed, each the median of 50 frames drawn at random from the middle 80% of the recording; every frame is scored against these candidates and the template refined from the top 10% by agreement, so movement-corrupted frames are suppressed. The result is cached next to the raw file.

`buildFusiAnatomy` then:

1. **Upsamples** every session once to isotropic atlas resolution and equalizes it to a common median, so no session dominates the average.
2. **Asks for a seed session** via `selectSeedSession` — a GUI showing every session as a column of coronal slices. The choice is cached in `seed_session_for_anatomy.txt`; on later runs the GUI is skipped.
3. **Rigidly registers** every other session to the seed with elastix, once, caching each transform under `anatomy_registration/<session>/`.
4. **Averages** the aligned sessions (NaN-omitting) into `<mouse>_anatomy.mat`.

The template stays in the **native session orientation** — only the sampling changes — so it feeds straight into `setupFusiOptions`.

```matlab
anatomy = buildFusiAnatomy(opts, sessionvols, pxsizesession, sessionnames);
```

| Output field | Contents |
| :--- | :--- |
| `volume` | `[D1 D2 D3]` anatomy at atlas resolution |
| `voxelsize_mm` | resolution of `volume` (`== opts.atlas_res`) |
| `seed_session` | index of the seed |
| `tforms{i}` | session *i* → seed rigid transform — **reused in step 5**, never re-fitted |
| `sessionvolumes` | every session warped into seed space (QC) |
| `session_names`, `session_voxelsize_mm` | bookkeeping |

**Choosing a seed matters.** It defines the common space that everything downstream is expressed in, so pick a session that is well centred, artefact-free, and covers as much of the brain as possible. A bad seed costs you at every later step.

**Averaging vs. the seed alone.** The mean is the robust default but blurs vessels wherever the alignment is imperfect; the seed session on its own is sharper. Both are valid registration targets — `anatomy.sessionvolumes(:,:,:,anatomy.seed_session)` gives you the second option.

**A harder alternative.** `buildFusiAnatomyNaive` runs the same idea *groupwise* over several passes: the template is rebuilt each pass as the median of the sessions that passed an alignment check (enough overlap, positive correlation), and each fit is warm-started from the previous pass. Use it when one or two sessions are bad enough to poison a single-pass average — it leaves them out instead of aborting.

---

## 2. Brain Orientation

**Function:** `setupFusiOptions(volume, voxelsize_mm, opts)`

Tells LightSuite how the sample axes map onto the atlas axes. On the first run a GUI appears; the chosen permutation is written to `brain_orientation.txt` and reused silently afterwards. The function also attaches the anatomy to `opts` as the registration sample and writes `regopts.mat`, which every later step reads.

```matlab
opts = setupFusiOptions(anatomy.volume, opts.atlas_res, opts);
```

> **Pass the atlas resolution, not the raw one.** The anatomy from step 1 is **already at atlas resolution**. Passing `pxsizesession` here rescales it a second time, and every subsequent step is wrong by that factor.

---

## 3. Place Matched Control Points

**Function:** `matchControlPoints_minimal(opts)`

A side-by-side GUI — your anatomy on the left, the vascular atlas on the right, both sliced along dimension 1. Click a landmark on one side, then its counterpart on the other. Vessel bifurcations, the midline, the brain outline and the ventricles all make good landmarks.

Each sample slice is presented **four times**, each showing about 60% of the image with one side blacked out (left, right, top, bottom). Every presented half is an independent data point: it gets its own atlas slice, taken in order from its own points, from interpolation between its own anchors, from the running affine fit (once ≥5 non-coplanar pairs exist), or from the nearest anchor's offset. Halves never borrow each other's slice, so a saved session always reopens on exactly the atlas slice it was annotated on.

### Controls

| Key / Action | Function |
| :--- | :--- |
| **← / →** | Previous / next presented half |
| **Scroll** | Move the atlas slice |
| **Click** | Add a point (alternate sample and atlas) |
| **Backspace** | Delete the most recently added point |
| **c** | Clear both point sets on the current half |
| **Space** | Toggle the annotation overlay |
| **Enter** | Jump to a sample slice number |
| **1–9 / 0** | Sample channel (0 = RGB from the first three) |
| **q** | Set the sample display quantiles (contrast) |
| **s** | Save to `control_points_minimal.mat` |

Points are reloaded when you reopen the GUI on the same folder, so annotation can be spread over several sittings.

**How many?** At least five pairs spanning all three dimensions are required. The paper used a median of 219 per mouse; spread over the full anteroposterior extent matters more than raw count, since a dense cluster constrains the fit only where it sits.

---

## 4. Fit the Anatomy → Atlas Transform

**Function:** `multiobjRegistrationFusi(opts, contol_point_wt, usemultistep)`

Fits an affine and then a B-spline transform, optimizing a weighted sum of **image similarity** (mutual information between the vascular atlas and your anatomy, after a high-pass that isolates vessel contrast) and the **distance between your matched control points**. Both stages see the landmark term, so it steers the whole chain rather than correcting it afterwards.

```matlab
wtpoints                   = 0.1;   % landmark weight vs image similarity
opts.bspline_spatial_scale = 1.6;   % mm — smaller = more local deformation
opts.n_histogram_bins      = 48;    % bins for the mutual-information estimate
multiobjRegistrationFusi(opts, wtpoints, false);
```

| Parameter | Meaning | Tuning |
| :--- | :--- | :--- |
| `contol_pt_wt` | Weight of the landmark term | Raise (`0.1 → 1 → 2`) if the fit drifts away from your points |
| `bspline_spatial_scale` | B-spline control-point spacing, mm | Raise if the result looks over-warped; lower for more local deformation |
| `n_histogram_bins` | Bins for mutual information | 48 is a good default; lower for very noisy anatomies |
| `usemultistep` | Multi-resolution B-spline schedule | `false` for most fUS data |

The function writes per-dimension overlay PNGs (`<mouse>_dim*_affine_registration.png`, `<mouse>_dim*_bspline_registration.png`) at both stages. **Look at them before trusting anything downstream** — a bad registration is not detectable from the maps themselves.

`transform_params.mat` holds the resulting chain:

| Field | Contents |
| :--- | :--- |
| `tform_affine_samp20um_to_atlas_10um_px` | Affine sample → atlas |
| `tform_bspline_samp20um_to_atlas_20um_px` | Path to the elastix B-spline parameter file |
| `atlassize`, `atlasres` | Atlas grid the warps target |
| `ori_size`, `ori_pxsize`, `how_to_perm` | The sample the chain was fit on |

`transformFusiAnnotationVolume(opts)` runs the chain **backwards**, bringing the Allen annotation into your anatomy's own space — the quickest way to see which areas your slab actually covered.

---

## 5. Applying the Transform to Functional Data

### 5a. Recording → anatomy (rigid)

**Function:** `fusiRecordingToAnatomy(opts, data, voxelsize_mm, recname)`

A functional recording sits wherever the probe was on that day, so it needs its own rigid step into seed space. This fits it against the recording's time-median, caches the transform under `<savepath>/<recname>/rigid_atlas_to_samp_20um.txt`, and prints a QC overlay.

```matlab
tformfuntoanatomy = fusiRecordingToAnatomy(opts, scandata, pxsizesession, recname);
```

If the recording is **one of the sessions that built the anatomy**, its rigid transform is already known and re-fitting it only adds noise — reuse it directly:

```matlab
tformfuntoanatomy = struct('tformrigid', anatomy.tforms{irec}, ...
                           'Ranatomy',   imref3d(size(anatomy.volume)));
```

### 5b. Anything → atlas space

**Function:** `applyFusiTransforms(opts, data, voxelsize_mm, savefullvols, transfuntoanatomy)`

Chains the rigid step with the affine + B-spline warp and returns the result in atlas space. `data` may be a single 3-D map or a 4-D `[d1 d2 d3 nframes]` array; the dense B-spline displacement field is computed **once** and reused for every frame, so a whole recording costs one transformix call plus one `imwarp` per frame.

```matlab
res = applyFusiTransforms(opts, outmap.Data, pxsizesession, true, tformfuntoanatomy);
```

| Output field | Contents |
| :--- | :--- |
| `vreg` | `[atlassize × nframes]` — full non-rigid warp (only if `savefullvols`) |
| `vregaff` | the same, **affine only** — a free sanity check on the B-spline |
| `areasignals` | `[nArea × nframes]` mean signal per parcellation leaf, non-rigid |
| `areasignalsaff` | the same, affine only |
| `groupidx` | Allen parcellation index labelling the rows above |
| `areavols` | volume of each area, in mm³ |

> **Set `savefullvols = false` unless you need the voxels.** Full atlas-space movies are `264 × 160 × 228 × nframes` singles. The per-area timecourses in `areasignals` are what most analyses use, and they cost nothing to keep.

**Reading `areasignals` by area name.** Rows are parcellation *leaves*; a named structure such as `VISp` spans several of them (its layers). Select the level you want in `parcelinfo`, then pool the matching rows weighted by volume:

```matlab
isstruct = strcmp(opts.parcelinfo.parcellation_term_set_name, 'structure');
leafidx  = opts.parcelinfo.parcellation_index( ...
    isstruct & strcmpi(opts.parcelinfo.parcellation_term_acronym, 'VISp'));
irows    = ismember(res.groupidx, leafidx);
tc       = sum(res.areavols(irows) .* res.areasignals(irows, :), 1, 'omitnan') ...
           / sum(res.areavols(irows));
```

`plotFusiActivationMontage(vol, opts)` draws an atlas-space volume as a grid of coronal slices over a greyscale underlay (`'sliceaxis', 1` in CCF space, where AP runs along dimension 1).

---

## 6. Cortical Flatmaps

**Function:** `fusiVolumeToFlatmap(vol, opts)`

Any volume already in CCF space can be flattened onto the Allen cortical surface: for every point of the flatmap, the volume is averaged (`'kind','mean'`) or maximized (`'max'`) along the corresponding cortical streamline.

MATLAB does **not** do the projection. It writes the volume to a small v7 `.mat`, calls `src/fusi/py/fusi_flatmap_project.py` with the interpreter in `opts.pythonexe`, and reads the 2-D result back. The streamlines live in `allensdk` / `ccf_streamlines`, so that environment can be any conda/venv — or another machine — without linking into MATLAB. The worker upsamples to the 10 µm grid the streamlines are defined on, so no resampling is needed on the MATLAB side.

```matlab
fmopts = struct();
fmopts.pythonexe   = 'C:\miniconda3\envs\allensdk\python.exe';
fmopts.resourcedir = 'C:\AllenAtlas\flatmaps';
fmopts.inputres_um = 50;        % atlas-space maps are 50 um isotropic
fmopts.hemisphere  = 'left';
fmopts.kind        = 'mean';

fm = fusiVolumeToFlatmap(res.vreg, fmopts);

imagesc(fm.flatmap.'); axis image off; colorbar
hold on; [yy, xx] = find(fm.boundaries.'); plot(xx, yy, '.k', 'MarkerSize', 1);
```

| Output field | Contents |
| :--- | :--- |
| `flatmap` | 2-D single, the streamline-reduced cortex |
| `boundaries` | area-boundary overlay image, same size (or `[]`) |
| `regionBoundaries` | `.names` (acronyms) and `.coords` (`N×2` `[row col]`), to outline a chosen subset of areas |
| `hemisphere`, `kind`, `inputres_um`, `maskres_um` | the settings used |

**Symmetrizing first.** fUS coverage is rarely identical on the two sides of the brain, and the flatmap is drawn per hemisphere. Mirror-averaging across the midline before projecting fills one hemisphere with whatever the other one saw:

```matlab
Nmid = size(vol, 3)/2;
vol  = (vol(:, :, 1:Nmid) + flip(vol(:, :, Nmid+1:end), 3))/2;
vol  = cat(3, vol, flip(vol, 3));
```

Skip it whenever left/right differences are part of the result. `symmetrizeVol` / `symmetrizeVolNan` do the same thing, NaN-aware.


## Output File Summary

```
<mouse>/lightsuite/
├── seed_session_for_anatomy.txt     # cached seed choice
├── anatomy_registration/<session>/  # cached session -> seed rigid transforms
├── <mouse>_anatomy.mat              # the anatomy scan + all session transforms
├── brain_orientation.txt            # cached sample -> atlas axis permutation
├── regopts.mat                      # atlas, annotation, sample, resolutions
├── control_points_minimal.mat       # your matched landmarks
├── transform_params.mat             # the fitted affine + B-spline chain
├── bspline_samp_to_atlas_20um.txt   # inverted B-spline (atlas -> sample)
├── <mouse>_dim*_affine_registration.png
├── <mouse>_dim*_bspline_registration.png
└── <recname>/
    ├── rigid_atlas_to_samp_20um.txt # recording -> anatomy rigid transform
    ├── dim1_rigid_registration.png
    ├── corr_registered.mat          # a map in atlas space
    └── time_registered.mat          # per-area timecourses in atlas space
```

---

## Troubleshooting

| Symptom | Likely cause |
| :--- | :--- |
| Anatomy looks like a blurred smear | Sessions did not align to the seed. Check `anatomy.sessionvolumes` slice by slice; try `buildFusiAnatomyNaive`, or a better seed (delete `seed_session_for_anatomy.txt` to re-pick). |
| Registration overlays are grossly off | Wrong axis permutation. Delete `brain_orientation.txt` and re-run `setupFusiOptions`. |
| Fit ignores your landmarks | Raise `contol_pt_wt`. Also check the GUI reported the pairs — halves with unequal sample/atlas point counts are silently dropped by `multiobjRegistrationFusi`. |
| Result looks over-warped | Raise `bspline_spatial_scale`; compare `vreg` against `vregaff` to see how much the B-spline added. |
| Out of memory in `applyFusiTransforms` | Pass `savefullvols = false` and use `areasignals`. |
| Flatmap call fails | `opts.pythonexe` does not point at an environment with `allensdk` + `ccf_streamlines`, or `resourcedir` is missing one of the four resource files. |

---

## References

* Atlas coordinate framework: Allen Mouse Brain Common Coordinate Framework (Wang et al., *Cell*, 2020).
* Cortical flatmaps: [ccf_streamlines](https://github.com/AllenInstitute/ccf_streamlines) (Allen Institute).
* `mapCorrelation` and `hemodynamicResponse` follow the fUS analysis conventions of the Urban and Macé labs (Montaldo, Macé et al.).
* LightSuite: Karamanlis et al., *Scalable atlas registration and cell detection across the central nervous system* (in preparation).
