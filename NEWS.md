# rgeoadaptels 0.2.0

# rgeoadaptels 0.3.1

## Fixed

* **An `sf` point layer passed to `grow_seeds_raster()` failed with a confusing
  message about column counts.** `sf` and `sfc` inherit from `data.frame`, so
  they matched the matrix branch before the spatial one, and their geometry
  list-column was read as coordinates -- "got 1 column(s)". The spatial classes
  are now tested first.

  The ordering is not cosmetic. The matrix branch has no CRS, so it takes
  coordinates on faith as already being in the raster's system; the spatial
  branch reprojects and says so. A frame that genuinely held x/y columns in
  another system would have been read silently wrong, which is the worst kind
  of wrong. `?grow_seeds_raster` now states that plainly.


# rgeoadaptels 0.3.0

## Added

* **`grow_seeds` / `grow_seeds_raster`: seeded spectral region growing
  ("inverse OBIA").** Each operator-placed point is grown into the region that
  looks like the pixel it sits on; everything unseeded is left unassigned
  (`-1`). The inverse of `adaptels()`/`sicle()`, built to delineate standing
  dead trees from a hand-digitised point layer.

  It calls the SICLE IFT kernel (`C_ift_fmax`) once, with every seed, and
  never removes one, so `labels == i` is the region grown from the i-th point.
  Options: `max_cost` (a ΔE tolerance on CIELAB input), `band_weights`,
  `compactness`, `seed_window`, `max_radius`, `fill_holes`. The kernel is
  unchanged, so this is the R twin of plGeoAdaptels' `grow_seeds`; the two are
  bit-identical on 30 cases in `tools/cross_validate_against_plgeoadaptels.R`.
  The cleanup helpers (`fill_holes`) are pure R, matching scipy on the Python
  side. `grow_seeds_raster` reads any terra-readable point layer, reprojects
  it when needed, and writes a label raster plus crown polygons. See
  `docs/grow_seeds_guide.md`.

* The cross-check now pins plGeoAdaptels to `v0.9.0`, the release that adds
  `grow_seeds`.


* **`sicle()`** -- Superpixels through Iterative CLEarcutting. Oversample
  seeds, grow an optimum-path forest with the `fmax` path-cost and the
  `wroot` arc-cost, score every seed, discard the least relevant, repeat.
  Optional saliency map, optional mask, explicit or sampled seeds.
* Still bit-identical to plGeoAdaptels, now across **20** cross-validation
  cases: the ten adaptel cases, three `enforce_connectivity` cases, and
  seven SICLE cases covering three iteration counts, three segment counts,
  a saliency map and a mask.

## Two things worth knowing about SICLE

**Seeds are the one part that does not match by default.** With
`seeds = NULL` they are drawn with R's own RNG, so `set.seed()` controls
them, and they will not be the seeds plGeoAdaptels draws: NumPy's
`Generator.choice` cannot be reproduced outside NumPy, and reimplementing
an undocumented ordering detail of a third-party library is what left rHRG
disagreeing with `scikit-image`'s watershed. Belem et al. treat the
sampling as a free choice rather than part of the algorithm, so passing
`seeds` takes it out of the comparison -- which is how the cross-check gets
to be an equality.

**`n_iterations` does less than it looks.** The paper's curve is
`M(i) = max(N0^(1 - i/(Omega-1)), Nf)`, whose exponent never sees the
target, so `n_iterations = 3` can be bit-identical to `2`. That is the
paper's design, not a defect here, and there is a test pinning it.

# rgeoadaptels 0.1.0

First release. The R implementation of plGeoAdaptels.

* `adaptels()` grows scale-adaptive superpixels from raster bands. Minkowski,
  cosine and angular distances, 4- or 8-connectivity, an optional per-band
  normalisation, and a nodata mask.
* `enforce_connectivity()` splits any adaptel that arrived in more than one
  piece, optionally absorbing slivers below `min_size`.
* `read_bands()` and `adaptels_raster()` bridge to \pkg{terra}, which is only
  a Suggests: the segmentation works on plain arrays with no spatial
  dependency.
* R and plain C, no Rcpp, no required packages.
* **Bit-identical to plGeoAdaptels 0.5.0** across all thirteen cases in
  `tools/cross_validate_against_plgeoadaptels.R` — every metric, both
  connectivities, a non-default Minkowski exponent, the normalise path, a
  mask, single- and multi-band input, and `enforce_connectivity` at three
  `min_size` values. Also checked locally against the real 400x400 scene in
  the plGeoAdaptels repository, which is not shipped here.

  Unlike the other twin pairs in this family, the cross-check is an equality
  and not a tolerance. rHRG reimplements scikit-image's watershed and differs
  where plateau ties fall; GeoPaletteR computes in double where GeoPalette
  stores single precision. Neither applies to a port of the same kernel.
