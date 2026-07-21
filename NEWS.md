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
