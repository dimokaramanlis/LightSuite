#!/usr/bin/env python
"""Project a fUSI atlas-space volume onto an Allen cortical flatmap.

Standalone worker called by the MATLAB function fusiVolumeToFlatmap.m (via
system() / the env's python.exe). It mirrors the projection used in
generate_flatmaps.ipynb but for whole-cortex (not per-layer) averaging of a
single arbitrary volume:

  1) read the input volume (CCF orientation AP x DV x ML, isotropic) from a
     v7 .mat file;
  2) upsample it to the 10 um CCF grid (1320, 800, 1140) that the
     ccf_streamlines surface paths are defined on;
  3) build a whole-isocortex mask from the Allen annotation (fetched with the
     allensdk at a chosen resolution, default 50 um to match the data, then
     nearest-upsampled to 10 um and cached);
  4) project along cortical streamlines with ccf_streamlines
     (Isocortex2dProjector.project_volume, kind='mean' by default);
  5) write the 2-D flatmap (+ optional area-boundary image) to an output .mat.

Everything is exchanged as v7 .mat so scipy.io can read/write it (scipy cannot
read MATLAB v7.3/HDF5), and the volume is passed at its native (small)
resolution so only a ~40 MB file crosses between MATLAB and Python.

Usage:
  python fusi_flatmap_project.py --input in.mat --output out.mat \
      --resourcedir D:/AllenAtlas/flatmaps [--hemisphere left] [--kind mean] \
      [--maskres 50] [--no-boundaries]
"""
import argparse
import os
import sys

import numpy as np
from scipy.io import loadmat, savemat
from scipy.ndimage import zoom

# 10 um CCF full-brain grid the streamlines / surface paths live on.
TARGET_SHAPE = (1320, 800, 1140)   # (AP, DV, ML)
ISOCORTEX_ID = 315                 # Allen structure id for Isocortex


# --------------------------------------------------------------------------
def _fit_shape(vol, shape):
    """Crop or zero/false-pad `vol` to exactly `shape` (guards zoom rounding)."""
    if vol.shape == tuple(shape):
        return vol
    out = np.zeros(shape, dtype=vol.dtype)
    sl = tuple(slice(0, min(a, b)) for a, b in zip(vol.shape, shape))
    out[sl] = vol[sl]
    return out


def _upsample(vol, shape, order):
    """Resample `vol` onto `shape`. Integer factors use exact np.repeat."""
    if vol.shape == tuple(shape):
        return vol
    factors = [t / s for t, s in zip(shape, vol.shape)]
    if all(float(f).is_integer() for f in factors):
        r = [int(round(f)) for f in factors]
        out = np.repeat(np.repeat(np.repeat(vol, r[0], 0), r[1], 1), r[2], 2)
        return _fit_shape(out, shape)
    out = zoom(vol, factors, order=order)
    return _fit_shape(out, shape)


def load_or_build_cortex_mask(resourcedir, maskres, target_shape):
    """Whole-isocortex boolean mask on the 10 um grid.

    The annotation is fetched/built at `maskres` um (small, matching the fUSI
    resolution), cached there, then nearest-upsampled to the 10 um grid.
    """
    cache = os.path.join(resourcedir, "cortex_mask_%dum.npy" % int(maskres))
    if os.path.exists(cache):
        maskn = np.load(cache)
    else:
        from allensdk.core.reference_space_cache import ReferenceSpaceCache
        print("  building isocortex mask from allensdk (%d um)..." % int(maskres))
        rspc = ReferenceSpaceCache(
            resolution=int(maskres),
            reference_space_key="annotation/ccf_2017",
            manifest=os.path.join(resourcedir, "manifest.json"),
        )
        annotation, _ = rspc.get_annotation_volume()
        tree = rspc.get_structure_tree(structure_graph_id=1)
        ids = tree.descendant_ids([ISOCORTEX_ID])[0]
        maskn = np.isin(annotation, ids)
        del annotation
        np.save(cache, maskn)
    return _upsample(maskn, target_shape, order=0).astype(bool)


def compute_region_boundaries(resourcedir):
    """Allen area boundaries in flatmap pixels, keyed by region (dict acronym->Nx2)."""
    import ccf_streamlines.projection as ccfproj
    bf = ccfproj.BoundaryFinder(
        projected_atlas_file=os.path.join(resourcedir, "flatmap_butterfly.nrrd"),
        labels_file=os.path.join(resourcedir, "labelDescription_ITKSNAPColor.txt"),
    )
    return bf.region_boundaries()


def boundaries_to_image(boundaries, flat_shape):
    """Rasterize a region-boundary dict onto a flat_shape uint8 image."""
    img = np.zeros(flat_shape, dtype=np.uint8)
    for _, coords in boundaries.items():
        c = np.round(coords).astype(int)
        valid = ((c[:, 0] >= 0) & (c[:, 0] < flat_shape[0]) &
                 (c[:, 1] >= 0) & (c[:, 1] < flat_shape[1]))
        c = c[valid]
        img[c[:, 0], c[:, 1]] = 1
    return img


# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="input v7 .mat (vol, inputres_um)")
    p.add_argument("--output", required=True, help="output .mat to write")
    p.add_argument("--resourcedir", required=True,
                   help="dir with flatmap_butterfly.h5/.nrrd, surface_paths_10_v3.h5, "
                        "labelDescription_ITKSNAPColor.txt, manifest.json")
    p.add_argument("--hemisphere", default="left", choices=["left", "right"])
    p.add_argument("--kind", default="mean", choices=["mean", "max"])
    p.add_argument("--maskres", type=int, default=50,
                   help="annotation resolution (um) for the cortex mask (default 50)")
    p.add_argument("--no-boundaries", action="store_true",
                   help="skip the area-boundary overlay image")
    args = p.parse_args()

    import ccf_streamlines.projection as ccfproj

    # ---- input volume -----------------------------------------------------
    md = loadmat(args.input)
    if "vol" not in md:
        sys.exit("input .mat has no 'vol' variable")
    vol = np.ascontiguousarray(np.squeeze(md["vol"])).astype(np.float32)
    inputres = float(np.squeeze(md["inputres_um"])) if "inputres_um" in md else 50.0
    print("  input volume %s at %g um -> 10 um grid %s" %
          (vol.shape, inputres, TARGET_SHAPE))
    if vol.ndim != 3:
        sys.exit("expected a 3-D volume, got shape %s" % (vol.shape,))

    # ---- upsample data to the streamline (10 um) grid --------------------
    vol10 = _upsample(vol, TARGET_SHAPE, order=1)
    del vol

    # ---- whole-cortex mask ------------------------------------------------
    mask = load_or_build_cortex_mask(args.resourcedir, args.maskres, TARGET_SHAPE)
    if mask.shape != vol10.shape:
        sys.exit("mask shape %s != volume shape %s" % (mask.shape, vol10.shape))

    # neutralise non-cortex voxels in place (no second full-size buffer):
    #   mean -> NaN drops out of nanmean; max -> float32 min never wins.
    if args.kind == "mean":
        vol10[~mask] = np.nan
    else:
        vol10[~mask] = np.finfo(np.float32).min
    del mask

    # ---- project along cortical streamlines ------------------------------
    proj = ccfproj.Isocortex2dProjector(
        os.path.join(args.resourcedir, "flatmap_butterfly.h5"),
        os.path.join(args.resourcedir, "surface_paths_10_v3.h5"),
        hemisphere=args.hemisphere,   # single hemisphere (no butterfly view)
    )
    print("  projecting whole cortex (%s, kind=%s)..." % (args.hemisphere, args.kind))
    with np.errstate(all="ignore"):
        flat = proj.project_volume(vol10, kind=args.kind)
    flat = np.asarray(flat, dtype=np.float32)
    if args.kind == "mean":
        flat = np.nan_to_num(flat)
    else:
        f32min = np.finfo(np.float32).min
        flat = np.where(np.isfinite(flat) & (flat > f32min), flat, 0.0).astype(np.float32)
    del vol10

    out = {"flatmap": flat,
           "flat_shape": np.asarray(flat.shape, dtype=np.float64),
           "hemisphere": args.hemisphere,
           "kind": args.kind,
           "inputres_um": np.float64(inputres),
           "maskres_um": np.float64(args.maskres)}

    # ---- area boundaries (optional overlay) ------------------------------
    # `boundaries`      = merged uint8 outline image (all areas, back-compat)
    # `region_names`    = 1xN object array of region acronyms/labels
    # `region_coords`   = 1xN object array of Nx2 [row col] flatmap-pixel arrays
    # so MATLAB can colour a chosen subset of areas (e.g. the dorsal/ventral
    # visual streams) rather than only the merged overlay.
    if not args.no_boundaries:
        try:
            bnd = compute_region_boundaries(args.resourcedir)
            out["boundaries"] = boundaries_to_image(bnd, flat.shape)
            names = list(bnd.keys())
            coords = [np.asarray(bnd[k], dtype=np.float64) for k in names]
            out["region_names"] = np.array(names, dtype=object)
            out["region_coords"] = np.array(coords + [None], dtype=object)[:-1]
        except Exception as exc:                       # non-fatal overlay
            print("  (boundary image skipped: %s)" % exc)
            out["boundaries"] = np.zeros(flat.shape, dtype=np.uint8)

    savemat(args.output, out, do_compression=True)
    print("  wrote %s (flatmap %s)" % (args.output, flat.shape))


if __name__ == "__main__":
    main()
