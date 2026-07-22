# rgeoadaptels

<img src="https://raw.githubusercontent.com/igorpawelec/rgeoadaptels/main/www/rgeoadaptels.png" alt="rgeoadaptels logo" align="right" width="200"/>

[![R-CMD-check](https://github.com/igorpawelec/rgeoadaptels/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igorpawelec/rgeoadaptels/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%203.6-blue.svg)](https://www.r-project.org)

**Scale-adaptive superpixels for geospatial raster data.**

R and plain C. No Rcpp, no required packages, and `terra` only if you want to
read files.

> **Python users:** the same algorithm is in
> [plGeoAdaptels](https://github.com/igorpawelec/plGeoAdaptels). The two are
> separate repositories because their tooling and idioms do not mix, but they
> are **bit-identical** — see [Agreement](#agreement-with-plgeoadaptels).

## Install

```r
# install.packages("remotes")
remotes::install_github("igorpawelec/rgeoadaptels")
```

## Use

```r
library(rgeoadaptels)

# data is (bands, rows, cols), or a plain matrix for one band
seg <- adaptels(data, threshold = 60)
seg$n_adaptels
seg$labels          # integer matrix, 0-based ids, -9999 for nodata
```

There is no target count. Regions grow from a seed until their internal
distance passes `threshold`, and the pixels beyond it become seeds in turn,
so the number of adaptels follows the scene.

With a raster:

```r
adaptels_raster("scene.tif", "adaptels.tif", threshold = 60)
```

### The threshold is per metric, not universal

This is the single easiest way to get a nonsense result, so the package
refuses rather than obliges:

| distance | range | sensible threshold |
|---|---|---|
| `minkowski` | grows with the data range | 60 on 0-255 imagery |
| `cosine` | [0, 1] by construction | around 0.03 |
| `angular` | [0, 1] by construction | around 0.03 |

Passing 60 to `cosine` would merge the whole raster into one adaptel. It
raises instead, and says what to try. The same applies to `normalize = TRUE`:
normalising caps the largest possible Minkowski distance at
`n_bands^(1/p)` — about 1.73 for three bands — so the raw-data default cannot
be reached and is rejected.

### Adaptels can arrive in more than one piece

Adaptels compete: a later one takes a pixel from an earlier one whenever it
arrives with a smaller accumulated distance. That competition is what gives
the method its boundary adherence, and it can also cut an earlier adaptel in
two. On a 400x400 three-band scene at the default threshold, about 10 per
cent come out split.

Harmless if the labels are a lookup. Not harmless for zonal statistics, which
would average two spatially separate patches into one "object":

```r
seg   <- adaptels(data, threshold = 60)
split <- enforce_connectivity(seg$labels)     # every region now contiguous
```

Not applied automatically, because it changes the adaptel count.

## Agreement with plGeoAdaptels

Bit-identical, and checked rather than asserted:

```sh
pip install plgeoadaptels
python3 tools/generate_plgeoadaptels_reference.py
Rscript tools/cross_validate_against_plgeoadaptels.R
```

Thirteen cases — every metric, both connectivities, a non-default Minkowski
exponent, the normalise path, a mask with an interior hole, single- and
multi-band input, a constant raster, and `enforce_connectivity` at three
`min_size` values. Zero differing pixels in all of them. The real 400x400
three-band scene in the plGeoAdaptels repository was checked the same way
locally; it is not shipped here.

The check is an **equality**, not a tolerance, and that is worth a note
because the other twin pairs in this family could not manage it. rHRG
reimplements `scikit-image`'s watershed and differs on 0.25 % of pixels where
plateau ties fall differently. GeoPaletteR computes in double where
GeoPalette stores single precision. Neither problem exists here: this is a
port of the same kernel — same heap, same neighbour order, same arithmetic —
so there is no second implementation to disagree with, and anything short of
identical would be a bug.

Three details make that possible, and all three are easy to get wrong:

- **R stores matrices by column, Python by row.** The C is handed
  `t(matrix)`, whose column-major buffer is byte-identical to a row-major
  NumPy array. Skip it and the neighbour order transposes, which changes
  which adaptel claims a contested pixel.
- **Bands are indexed `layers[l * size + i]`**, matching NumPy's
  `(n_layers, size)` C-order rather than an R matrix's layout.
- **The heap is 1-based**, as in the Python, so the sift arithmetic is the
  same expression rather than a translated one.

## Citation

Pawelec, I. (2026). *rgeoadaptels: Scale-Adaptive Superpixels for Geospatial
Raster Data*. R package version 0.1.0.
https://github.com/igorpawelec/rgeoadaptels

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

## Licence

GPL-3.
