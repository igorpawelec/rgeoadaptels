"""Write reference segmentations from plGeoAdaptels.

    pip install plgeoadaptels
    python3 tools/generate_plgeoadaptels_reference.py

Synthetic scenes only, so this runs on a machine that has no raster data.
The agreement was also checked locally against the real 400x400 three-band
scene in the plGeoAdaptels repository, which cannot be shipped here.

Covers every parameter path, because the ones left out are exactly the ones
that drift: connectivity 4 and 8, all three metrics on their own scales, a
non-default Minkowski exponent, the normalise path, single-band and
multi-band input, a mask, and enforce_connectivity at two min_size values.

Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
"""
import os

import numpy as np
from plgeoadaptels import adaptels_from_array, enforce_connectivity
from plgeoadaptels.sicle import sicle_from_array

OUT = os.path.join("tools", "reference")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(17)

# Deliberately non-square, so a transposed buffer is an error rather than a
# rearrangement that still fits.
SCENES = {
    "multi": rng.integers(0, 256, (3, 35, 48)).astype(np.float64),
    "single": rng.integers(0, 256, (1, 29, 33)).astype(np.float64),
    "flat": np.full((2, 20, 26), 42.0),
}

# A mask with a hole in the middle and a nodata border, so nodata is not
# merely a corner case at the edge.
mask = np.zeros((35, 48), dtype=np.uint8)
mask[0, :] = 1
mask[:, -1] = 1
mask[12:18, 20:28] = 1

CASES = [
    ("mink_rook", "multi", dict(threshold=60.0)),
    ("mink_queen", "multi", dict(threshold=60.0, queen_topology=True)),
    ("mink_p3", "multi", dict(threshold=60.0, minkowski_p=3.0)),
    ("cosine", "multi", dict(threshold=0.03, distance="cosine")),
    ("angular", "multi", dict(threshold=0.03, distance="angular")),
    ("normalize", "multi", dict(threshold=0.4, normalize=True)),
    ("masked", "multi", dict(threshold=60.0, mask=mask)),
    ("single", "single", dict(threshold=60.0)),
    ("flat", "flat", dict(threshold=60.0)),
    ("tight", "multi", dict(threshold=8.0)),
]

lines = []
for tag, scene, kw in CASES:
    data = SCENES[scene]
    labels, n = adaptels_from_array(data, **kw)
    np.savetxt(os.path.join(OUT, "out_%s.csv" % tag), labels,
               delimiter=",", fmt="%d")
    lines.append("%s,%s,%d" % (tag, scene, n))
    print("  %-11s %-7s -> %d adaptels" % (tag, scene, n))

# enforce_connectivity, fed from a real segmentation rather than made-up ids.
#
# Grown with queen_topology on purpose. enforce_connectivity tests
# 4-connectivity while this grows over 8, so adaptels come out that are
# 8-connected but not 4-connected and the splitting path actually runs. The
# obvious choice — the default rook segmentation above — splits nothing on
# random noise, so it would have exercised only the pass-through and proved
# very little.
base, _ = adaptels_from_array(SCENES["multi"], threshold=200.0,
                              queen_topology=True)
np.savetxt(os.path.join(OUT, "conn_input.csv"), base, delimiter=",", fmt="%d")
for ms in (0, 5, 20):
    out, n = enforce_connectivity(base, min_size=ms)
    np.savetxt(os.path.join(OUT, "conn_%d.csv" % ms), out, delimiter=",",
               fmt="%d")
    lines.append("conn%d,conn,%d" % (ms, n))
    print("  conn min_size=%-2d        -> %d adaptels" % (ms, n))

# ── SICLE ────────────────────────────────────────────────────────────
#
# Seeds are written out and handed to both implementations. NumPy's
# Generator.choice cannot be reproduced outside NumPy, so the sampler stays
# out of the comparison and the algorithm stays in it. plGeoAdaptels grew a
# `seeds` argument in 0.6.0 for exactly this.
sicle_rng = np.random.default_rng(23)
sc = SCENES["multi"]
n_pix = sc.shape[1] * sc.shape[2]
flat = sicle_rng.choice(np.arange(n_pix), size=700, replace=False)
SEEDS = np.stack([flat // sc.shape[2], flat % sc.shape[2]], axis=1)
np.savetxt(os.path.join(OUT, "sicle_seeds.csv"), SEEDS, delimiter=",", fmt="%d")

# A saliency map with a NaN-free interior; nodata in saliency is rejected by
# both sides, and that rejection is tested in the unit suites rather than here.
sal = np.zeros(sc.shape[1:], dtype=np.float64)
sal[10:25, 15:35] = 1.0
sal[26:30, 5:12] = 0.5
np.savetxt(os.path.join(OUT, "sicle_saliency.csv"), sal, delimiter=",",
           fmt="%.17g")

SICLE_CASES = [
    ("s_it2", dict(n_segments=60, n_iterations=2)),
    ("s_it5", dict(n_segments=60, n_iterations=5)),
    ("s_it10", dict(n_segments=60, n_iterations=10)),
    ("s_many", dict(n_segments=200, n_iterations=2)),
    ("s_few", dict(n_segments=5, n_iterations=2)),
    ("s_sal", dict(n_segments=60, n_iterations=2, saliency=sal)),
    ("s_mask", dict(n_segments=60, n_iterations=2, mask=mask)),
]
for tag, kw in SICLE_CASES:
    kw = dict(kw)
    seeds = SEEDS
    if "mask" in kw:
        # Drop seeds that land on nodata; both sides reject them, and the
        # point here is the algorithm, not the validation.
        keep = kw["mask"][SEEDS[:, 0], SEEDS[:, 1]] == 0
        seeds = SEEDS[keep]
        np.savetxt(os.path.join(OUT, "sicle_seeds_masked.csv"), seeds,
                   delimiter=",", fmt="%d")
    labels, n = sicle_from_array(sc, seeds=seeds, quiet=True, **kw)
    np.savetxt(os.path.join(OUT, "out_%s.csv" % tag), labels, delimiter=",",
               fmt="%d")
    lines.append("%s,sicle,%d" % (tag, n))
    print("  %-11s sicle   -> %d superpixels" % (tag, n))

for name, arr in SCENES.items():
    for l in range(arr.shape[0]):
        np.savetxt(os.path.join(OUT, "%s_band%d.csv" % (name, l)), arr[l],
                   delimiter=",", fmt="%.17g")
np.savetxt(os.path.join(OUT, "mask.csv"), mask, delimiter=",", fmt="%d")

with open(os.path.join(OUT, "cases.csv"), "w") as fh:
    fh.write("\n".join(lines) + "\n")

import plgeoadaptels
print("\nplgeoadaptels %s, %d cases written to %s/"
      % (plgeoadaptels.__version__, len(lines), OUT))
