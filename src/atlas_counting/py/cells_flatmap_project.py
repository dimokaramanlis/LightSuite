#!/usr/bin/env python
"""Project cell coordinates (Allen CCF) onto the Allen cortical flatmap.

Standalone worker called by the MATLAB function cellsToFlatmap.m. Unlike the
volume projection (fusi_flatmap_project.py), every cell keeps its own
continuous position: nothing is binned or averaged along streamlines. The two
hemispheres stay separate, on the two wings of the butterfly flatmap.

Voxel table (built once, cached in <resourcedir>/cells_flatmap_cache)
  The ccf_streamlines surface paths are the cortical streamlines, each a list
  of 10 um voxels from pia to white matter. Streamlines that own flatmap
  pixels (flatmap_butterfly.h5 "view lookup") sit at the mean of their pixels;
  the ~10% that own none take the position of the nearest owned streamline
  (by pial voxel), as ccf_streamlines does for coordinates. For every left-
  hemisphere voxel crossed by a streamline the table stores one
  representative streamline (the one closest to the voxel's mean flatmap
  position): its flatmap position, the voxel's depth along it (arc length from
  pia) and its id. Right-hemisphere streamlines are mirror copies and are not
  stored. The table is built with bounded memory (hits are bucketed to disk by
  AP slab) and is read back memory-mapped, so a query touches only the pages
  it needs.

How a cell is mapped
  Right-hemisphere cells are mirrored onto the left table and placed back on
  the right wing afterwards. Streamlines cross ~60% of the cortical voxels in a
  patchy pattern, so a cell takes the trilinear interpolation of flatmap
  position and depth over its 8 surrounding voxels when the crossed ones carry
  at least --minweight of the trilinear weight (weights renormalised), and
  otherwise a Gaussian-weighted average over the crossed voxels of its 3x3x3
  neighbourhood. Where the contributing voxels lie far apart on the flatmap
  (a tear in the flattening) the highest-weight voxel is used. Cells with no
  crossed voxel within the 3x3x3 neighbourhood are left unmapped (NaN). Layer,
  reference (layer-normalised) depth and local layer edges come from the
  cortical_layers_10_v2.h5 extents of the highest-weight voxel's streamline and
  avg_layer_depths.json.

Coordinates in the input are CCF microns (AP, DV, ML); voxel i (0-based) spans
[10*i, 10*(i+1)) um - the ccf_streamlines convention.

Usage:
  python cells_flatmap_project.py --input in.mat --output out.mat \
      --resourcedir D:/AllenAtlas/flatmaps [--tearpx 10] [--minweight 0.5] \
      [--sigma 0.75] [--batch 100000] [--layermaps] [--no-boundaries]
"""
import argparse
import json
import os
import shutil
import time

import h5py
import numpy as np
from scipy.io import loadmat, savemat

VOLUME_SHAPE = (1320, 800, 1140)   # 10 um CCF (AP, DV, ML)
RES_UM = 10.0
LAYER_KEYS = ["Isocortex layer 1", "Isocortex layer 2/3", "Isocortex layer 4",
              "Isocortex layer 5", "Isocortex layer 6a", "Isocortex layer 6b"]
LAYER_NAMES = ["L1", "L2/3", "L4", "L5", "L6a", "L6b"]
BUTTERFLY_SPACE = 184              # ccf_streamlines preset for 'flatmap_butterfly'
CACHE_DIR = "cells_flatmap_cache"
CACHE_VERSION = 1
XY_SCALE = 16.0                    # flatmap px stored as uint16 in 1/16 px
DEPTH_SCALE = 10.0                 # depth stored as uint16 in 0.1 um


def log(msg):
    print("  " + msg, flush=True)


# ==========================================================================
# cache building (once)
# ==========================================================================
def build_streamline_table(resourcedir, cdir):
    """Left streamlines: path row, flatmap position, layer extents."""
    with h5py.File(os.path.join(resourcedir, "flatmap_butterfly.h5"), "r") as f:
        view_lookup = f["view lookup"][:]
        view_size = np.asarray(f.attrs["view size"][:], dtype=np.int64)

    surf = view_lookup[:, 1].astype(np.int64)
    order = np.argsort(surf)
    usurf, inv = np.unique(surf[order], return_inverse=True)
    prow = np.empty(usurf.size, dtype=np.int64)
    with h5py.File(os.path.join(resourcedir, "surface_paths_10_v3.h5"), "r") as f:
        # fancy indexing into this gzip dataset is ~1000x slower than reading
        # contiguous slabs, so read slabs and pick the needed entries
        d = f["volume lookup flat"]
        slab = 25_000_000
        for s0 in range(0, d.shape[0], slab):
            i0, i1 = np.searchsorted(usurf, [s0, s0 + slab])
            if i1 > i0:
                block = d[s0:min(s0 + slab, d.shape[0])]
                prow[i0:i1] = block[usurf[i0:i1] - s0]
                del block
        pial = f["paths"][:, 0].astype(np.int64)   # first voxel = pial voxel
    path_row = np.empty(surf.size, dtype=np.int64)
    path_row[order] = prow[inv]
    keep = path_row >= 0
    view_flat = view_lookup[keep, 0].astype(np.int64)
    path_row = path_row[keep]
    del view_lookup, surf, order, usurf, inv, prow

    pial_xyz = np.stack(np.unravel_index(pial, VOLUME_SHAPE), 1).astype(np.int32)
    del pial
    rows = np.flatnonzero(pial_xyz[:, 2] < VOLUME_SHAPE[2] // 2)
    rowpos = np.full(pial_xyz.shape[0], -1, np.int64)
    rowpos[rows] = np.arange(rows.size)

    px, py = np.unravel_index(view_flat, tuple(view_size))
    view_k = rowpos[path_row]
    n = np.bincount(view_k, minlength=rows.size)
    inview = n > 0
    fxy = np.zeros((rows.size, 2), np.float32)
    fxy[inview, 0] = np.bincount(view_k, weights=px, minlength=rows.size)[inview] / n[inview]
    fxy[inview, 1] = np.bincount(view_k, weights=py, minlength=rows.size)[inview] / n[inview]

    from scipy.spatial import cKDTree
    tree = cKDTree(pial_xyz[rows[inview]])
    _, nn = tree.query(pial_xyz[rows[~inview]])
    fxy[~inview] = fxy[inview][nn]
    del tree, pial_xyz
    log("%d left streamlines, %d own flatmap pixels" % (rows.size, inview.sum()))

    with h5py.File(os.path.join(resourcedir, "cortical_layers_10_v2.h5"), "r") as f:
        layer_ext = np.stack([f[k][:][rows] for k in LAYER_KEYS], 1).astype(np.float32)

    np.save(os.path.join(cdir, "stream_rows.npy"), rows)
    np.save(os.path.join(cdir, "stream_fxy.npy"), fxy)
    np.save(os.path.join(cdir, "stream_layer_ext.npy"), layer_ext)
    np.save(os.path.join(cdir, "view_flat.npy"), view_flat)
    np.save(os.path.join(cdir, "view_k.npy"), view_k.astype(np.int32))
    np.save(os.path.join(cdir, "view_size.npy"), view_size)
    return rows, fxy


def build_voxel_table(resourcedir, cdir, rows, fxy, rowchunk=20000, nbuckets=44):
    """Per crossed voxel: representative streamline, flatmap xy, depth."""
    tmp = os.path.join(cdir, "tmp_buckets")
    os.makedirs(tmp, exist_ok=True)
    rec = np.dtype([("vox", "<u4"), ("sid", "<i4"), ("depth", "<f4")])
    apw = int(np.ceil(VOLUME_SHAPE[0] / nbuckets))
    sid_of_row = np.full(int(rows.max()) + 1, -1, np.int32)
    sid_of_row[rows] = np.arange(rows.size, dtype=np.int32)
    nhits = 0
    fh = [open(os.path.join(tmp, "b%03d.bin" % b), "wb") for b in range(nbuckets)]
    try:
        with h5py.File(os.path.join(resourcedir, "surface_paths_10_v3.h5"), "r") as f:
            d = f["paths"]
            for r0 in range(0, min(d.shape[0], sid_of_row.size), rowchunk):
                r1 = min(r0 + rowchunk, d.shape[0], sid_of_row.size)
                sel = np.flatnonzero(sid_of_row[r0:r1] >= 0)
                if sel.size == 0:
                    continue
                P = d[r0:r1][sel]                                   # uint32
                valid = P > 0
                xyz = np.stack(np.unravel_index(P, VOLUME_SHAPE), -1).astype(np.int16)
                step = np.sqrt(np.sum(np.diff(xyz, axis=1).astype(np.float32) ** 2, 2))
                step[~(valid[:, 1:] & valid[:, :-1])] = 0
                arc = np.zeros(P.shape, np.float32)
                np.cumsum(step, 1, out=arc[:, 1:])
                arc *= RES_UM
                ii, kk = np.nonzero(valid)
                out = np.empty(ii.size, rec)
                out["vox"] = P[ii, kk]
                out["sid"] = sid_of_row[r0 + sel[ii]]
                out["depth"] = arc[ii, kk]
                b = xyz[ii, kk, 0] // apw
                order = np.argsort(b, kind="stable")
                cuts = np.searchsorted(b[order], np.arange(nbuckets + 1))
                for bb in range(nbuckets):
                    if cuts[bb + 1] > cuts[bb]:
                        out[order[cuts[bb]:cuts[bb + 1]]].tofile(fh[bb])
                nhits += ii.size
    finally:
        for h in fh:
            h.close()
    log("%d streamline voxels bucketed" % nhits)

    parts = []
    for bb in range(nbuckets):
        fn = os.path.join(tmp, "b%03d.bin" % bb)
        h = np.fromfile(fn, rec)
        os.remove(fn)
        if h.size == 0:
            continue
        h = h[np.argsort(h["vox"], kind="stable")]
        uv, start, cnt = np.unique(h["vox"], return_index=True, return_counts=True)
        g = np.repeat(np.arange(uv.size), cnt)
        fx = fxy[h["sid"], 0].astype(np.float64)
        fy = fxy[h["sid"], 1].astype(np.float64)
        mx = np.bincount(g, fx) / cnt
        my = np.bincount(g, fy) / cnt
        dist = (fx - mx[g]) ** 2 + (fy - my[g]) ** 2
        pick = np.lexsort((dist, g))[start]          # closest to the voxel mean
        ps = h["sid"][pick]
        pfile = os.path.join(tmp, "p%03d.npz" % bb)
        np.savez(pfile,
                 vox=uv.astype(np.uint32),
                 fx=np.round(fxy[ps, 0] * XY_SCALE).astype(np.uint16),
                 fy=np.round(fxy[ps, 1] * XY_SCALE).astype(np.uint16),
                 depth=np.round(np.minimum(h["depth"][pick] * DEPTH_SCALE, 65535)).astype(np.uint16),
                 sid=ps.astype(np.int32))
        parts.append((pfile, uv.size))
        del h, uv, start, cnt, g, fx, fy, mx, my, dist, pick, ps

    ntot = sum(n for _, n in parts)
    mm = {k: np.lib.format.open_memmap(os.path.join(cdir, "vox_%s.npy" % k), mode="w+",
                                       dtype=dt, shape=(ntot,))
          for k, dt in [("vox", np.uint32), ("fx", np.uint16), ("fy", np.uint16),
                        ("depth", np.uint16), ("sid", np.int32)]}
    i = 0
    for pfile, n in parts:        # buckets are AP-ordered -> voxel ids stay sorted
        with np.load(pfile) as z:
            for k in mm:
                mm[k][i:i + n] = z[k]
        os.remove(pfile)
        i += n
    for k in mm:
        mm[k].flush()
    del mm
    shutil.rmtree(tmp, ignore_errors=True)
    log("%d crossed voxels in the table" % ntot)


def ensure_cache(resourcedir):
    cdir = os.path.join(resourcedir, CACHE_DIR)
    vfile = os.path.join(cdir, "version.json")
    if os.path.exists(vfile):
        with open(vfile) as f:
            if json.load(f).get("version") == CACHE_VERSION:
                return cdir
    log("building voxel table in %s (first run only)..." % cdir)
    t0 = time.time()
    os.makedirs(cdir, exist_ok=True)
    rows, fxy = build_streamline_table(resourcedir, cdir)
    build_voxel_table(resourcedir, cdir, rows, fxy)
    with open(vfile, "w") as f:
        json.dump({"version": CACHE_VERSION}, f)
    log("voxel table built in %.0f s" % (time.time() - t0))
    return cdir


# ==========================================================================
# queries
# ==========================================================================
def reference_layers(resourcedir):
    """Reference (average) layer edges in um from pia: L1 top ... WM border."""
    with open(os.path.join(resourcedir, "avg_layer_depths.json")) as f:
        tops = json.load(f)
    return np.array([0.0, tops["2/3"], tops["4"], tops["5"], tops["6a"],
                     tops["6b"], tops["wm"]], dtype=np.float64)


def normalised_depth(depth, ext, ref_edges):
    """Depth along a streamline (um) -> layer index (1-6) and reference depth.

    ext: (n, 6, 3) [start end thickness] per layer; absent layers are all 0.
    Depths beyond the first / last present layer are clamped to it.
    """
    n = depth.shape[0]
    start, end = ext[:, :, 0], ext[:, :, 1]
    present = end > 0
    below = present & (depth[:, None] <= end)
    lay = np.where(below.any(1), np.argmax(below, 1),
                   ext.shape[1] - 1 - np.argmax(present[:, ::-1], 1))
    r = np.arange(n)
    s, e = start[r, lay], end[r, lay]
    frac = np.clip((depth - s) / np.maximum(e - s, 1e-6), 0.0, 1.0)
    nd = ref_edges[lay] + frac * (ref_edges[lay + 1] - ref_edges[lay])
    ok = present.any(1) & np.isfinite(depth)
    nd[~ok] = np.nan
    return np.where(ok, lay + 1, 0), nd


def map_batch(v, tab, args):
    """v: (n, 3) continuous left-hemisphere voxel coords (centres at integers)."""
    n = v.shape[0]
    base = np.round(v).astype(np.int32)
    offs = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1)
                     for k in (-1, 0, 1)], np.int32)
    nb = base[:, None, :] + offs[None]                          # (n, 27, 3)
    inb = np.all((nb >= 0) & (nb < np.array(VOLUME_SHAPE, np.int32)), 2)
    nb[~inb] = 0
    # voxel ids fit uint32; id 0 is path padding, so it never matches
    ind = ((nb[..., 0].astype(np.uint32) * np.uint32(VOLUME_SHAPE[1]) + nb[..., 1].astype(np.uint32))
           * np.uint32(VOLUME_SHAPE[2]) + nb[..., 2].astype(np.uint32))
    ind[~inb] = 0
    # one sorted lookup per batch (same dtype as the table: no array copies)
    uq, inv = np.unique(ind, return_inverse=True)
    vox = tab["vox"]
    pos = np.minimum(np.searchsorted(vox, uq), vox.size - 1)
    uhit = np.asarray(vox[pos]) == uq
    ph = pos[uhit]
    nu = uq.size
    ufx = np.zeros(nu, np.float32)
    ufy = np.zeros(nu, np.float32)
    udep = np.zeros(nu, np.float32)
    usid = np.full(nu, -1, np.int32)
    ufx[uhit] = tab["fx"][ph] / XY_SCALE
    ufy[uhit] = tab["fy"][ph] / XY_SCALE
    udep[uhit] = tab["depth"][ph] / DEPTH_SCALE
    usid[uhit] = tab["sid"][ph]
    inv = inv.reshape(ind.shape)
    hit = uhit[inv]
    fx, fy, dep, sid = ufx[inv], ufy[inv], udep[inv], usid[inv]

    dv = v[:, None, :].astype(np.float32) - nb
    w_tri = np.where(hit, np.prod(np.clip(1 - np.abs(dv), 0, None), 2), 0)
    w_gau = np.where(hit, np.exp(-np.sum(dv ** 2, 2) / (2 * args.sigma ** 2)), 0)
    del nb, dv, ind, inv

    fallback = w_tri.sum(1) < args.minweight
    w = np.where(fallback[:, None], w_gau, w_tri)
    wsum = w.sum(1)
    mapped = wsum > 0
    wn = w / np.maximum(wsum, 1e-12)[:, None]
    best = np.argmax(w, 1)
    r = np.arange(n)

    cx, cy, cd = (wn * fx).sum(1), (wn * fy).sum(1), (wn * dep).sum(1)
    used = w > 0
    spread = np.maximum(np.where(used, fx, -np.inf).max(1) - np.where(used, fx, np.inf).min(1),
                        np.where(used, fy, -np.inf).max(1) - np.where(used, fy, np.inf).min(1))
    tear = mapped & (spread > args.tearpx)
    cx[tear], cy[tear], cd[tear] = fx[r, best][tear], fy[r, best][tear], dep[r, best][tear]
    cx[~mapped] = cy[~mapped] = cd[~mapped] = np.nan
    return cx, cy, cd, sid[r, best], mapped, fallback & mapped, tear


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="input v7 .mat with coords_um (N x 3)")
    p.add_argument("--output", required=True)
    p.add_argument("--resourcedir", required=True)
    p.add_argument("--tearpx", type=float, default=10.0,
                   help="flatmap spread (px) of the contributing voxels above "
                        "which the highest-weight voxel is used")
    p.add_argument("--minweight", type=float, default=0.5,
                   help="trilinear weight the crossed voxels must carry; below "
                        "it the Gaussian 3x3x3 neighbourhood is used")
    p.add_argument("--sigma", type=float, default=0.75,
                   help="Gaussian width (voxels) of the 3x3x3 fallback")
    p.add_argument("--batch", type=int, default=100000, help="cells per batch")
    p.add_argument("--layermaps", action="store_true",
                   help="also return flatmaps of the local layer edges")
    p.add_argument("--no-boundaries", action="store_true")
    args = p.parse_args()
    t0 = time.time()

    coords = np.asarray(loadmat(args.input)["coords_um"], dtype=np.float64).reshape(-1, 3)
    ncells = coords.shape[0]
    log("%d cells" % ncells)

    cdir = ensure_cache(args.resourcedir)
    tab = {k: np.load(os.path.join(cdir, "vox_%s.npy" % k), mmap_mode="r")
           for k in ("vox", "fx", "fy", "depth", "sid")}
    layer_ext = np.load(os.path.join(cdir, "stream_layer_ext.npy"), mmap_mode="r")
    view_size = tuple(int(x) for x in np.load(os.path.join(cdir, "view_size.npy")))
    ref_edges = reference_layers(args.resourcedir)

    # ---- mirror right-hemisphere cells onto the left table ----------------
    mid_um = VOLUME_SHAPE[2] * RES_UM / 2
    finite = np.all(np.isfinite(coords), 1)
    isright = finite & (coords[:, 2] > mid_um)
    c = np.where(finite[:, None], coords, -1e6)
    c[isright, 2] = 2 * mid_um - c[isright, 2]
    v = c / RES_UM - 0.5                                   # voxel centres at integers

    flatxy = np.full((ncells, 2), np.nan)
    depth = np.full(ncells, np.nan)
    ndepth = np.full(ncells, np.nan)
    layer = np.zeros(ncells)
    edges = np.full((ncells, 7), np.nan)
    mapped = np.zeros(ncells, bool)
    fallback = np.zeros(ncells, bool)
    tear = np.zeros(ncells, bool)
    for b0 in range(0, ncells, args.batch):
        sl = slice(b0, min(b0 + args.batch, ncells))
        cx, cy, cd, bsid, m, fb, tr = map_batch(v[sl], tab, args)
        flatxy[sl, 0], flatxy[sl, 1], depth[sl] = cx, cy, cd
        mapped[sl], fallback[sl], tear[sl] = m, fb, tr
        idx = np.flatnonzero(m)
        if idx.size:
            usid, inv = np.unique(bsid[idx], return_inverse=True)
            ext = np.asarray(layer_ext[usid])[inv]            # (n, 6, 3)
            lay, nd = normalised_depth(cd[idx], ext, ref_edges)
            layer[b0 + idx], ndepth[b0 + idx] = lay, nd
            top = np.where(ext[:, :, 1] > 0, ext[:, :, 0], np.nan)
            edges[b0 + idx] = np.concatenate([top, ext[:, :, 1].max(1)[:, None]], 1)

    # ---- butterfly placement: right cells on the right wing ----------------
    # left wing: x in [0, maxx); right wing: x' = 2*maxx - x in (maxx, 2*maxx]
    maxx = view_size[0] - BUTTERFLY_SPACE
    flatxy[isright, 0] = 2 * maxx - flatxy[isright, 0]
    log("%d / %d cells mapped (%d trilinear, %d 3x3x3 fallback, %d at tears), %.1f s" %
        (mapped.sum(), ncells, (mapped & ~fallback).sum(), fallback.sum(), tear.sum(),
         time.time() - t0))

    out = {
        "flatxy": flatxy,                      # 0-based butterfly px [x y]
        "depth_um": depth,
        "depthnorm_um": ndepth,
        "layer": layer,
        "layeredges_um": edges,
        "mapped": mapped.astype(np.float64),
        "isright": isright.astype(np.float64),
        "fallback": fallback.astype(np.float64),
        "tear": tear.astype(np.float64),
        "ref_layeredges_um": ref_edges,
        "layer_names": np.array(LAYER_NAMES, dtype=object),
        "midline_x": np.float64(maxx),
    }

    # ---- cortex silhouette (and local layer edges) on the butterfly --------
    view_flat = np.load(os.path.join(cdir, "view_flat.npy"))
    px, py = np.unravel_index(view_flat, view_size)
    okx = px < maxx
    vm = np.zeros((2 * maxx + 1, view_size[1]), np.uint8)
    vm[px[okx], py[okx]] = 1
    vm[maxx + 1:] = vm[maxx - 1::-1]
    out["viewmask"] = vm
    if args.layermaps:
        view_k = np.load(os.path.join(cdir, "view_k.npy"))
        lm = np.full((2 * maxx + 1, view_size[1], 7), np.nan, np.float32)
        ext = np.asarray(layer_ext[view_k[okx]])
        lm[px[okx], py[okx]] = np.concatenate(
            [np.where(ext[:, :, 1] > 0, ext[:, :, 0], np.nan),
             ext[:, :, 1].max(1)[:, None]], 1)
        lm[maxx + 1:] = lm[maxx - 1::-1]
        out["layermaps_um"] = lm

    if not args.no_boundaries:
        try:
            import ccf_streamlines.projection as ccfproj
            bf = ccfproj.BoundaryFinder(
                projected_atlas_file=os.path.join(args.resourcedir, "flatmap_butterfly.nrrd"),
                labels_file=os.path.join(args.resourcedir, "labelDescription_ITKSNAPColor.txt"))
            names, coords_b, hemi = [], [], []
            for h, isr in (("left", 0.0), ("right_for_both", 1.0)):
                bnd = bf.region_boundaries(hemisphere=h,
                                           view_space_for_other_hemisphere="flatmap_butterfly")
                for k, val in bnd.items():
                    if np.asarray(val).size == 0:
                        continue
                    names.append(k)
                    coords_b.append(np.asarray(val, dtype=np.float64))
                    hemi.append(isr)
            out["region_names"] = np.array(names, dtype=object)
            out["region_coords"] = np.array(coords_b + [None], dtype=object)[:-1]
            out["region_isright"] = np.asarray(hemi)
        except Exception as exc:                      # non-fatal overlay
            log("(boundaries skipped: %s)" % exc)

    savemat(args.output, out, do_compression=True)
    log("wrote %s (%.1f s)" % (args.output, time.time() - t0))


if __name__ == "__main__":
    main()
