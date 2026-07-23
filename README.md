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

## SICLE superpixels

`adaptels` lets the scene decide the count; `sicle` lets you set it. It starts
from far more seeds than you want, grows an optimum-path forest, scores every
seed and discards the least relevant, and repeats until `n_segments` remain.

```r
seg <- sicle(data, n_segments = 200)
seg$labels            # 0-based ids, each a single 8-connected region
seg$n_superpixels
```

Seeds can be given instead of sampled. That is what lets the two twins be
compared on the algorithm rather than the sampler — `NumPy`'s `Generator.choice`
cannot be reproduced outside NumPy, so the sampling is kept out of the
comparison:

```r
sicle(data, seeds = cbind(rows, cols), n_segments = 200)   # 1-based (row, col)
```

`n_iterations` does less than it looks: the paper's preservation curve makes 3
bit-identical to 2. See `?sicle` for why, and why 2 is a speed setting rather
than a quality one.

## grow_seeds — seeded region growing

The inverse of `adaptels` and `sicle`. Rather than partition the whole image,
you supply points and each grows into the region that looks like the pixel it
sits on; everything unseeded stays `-1`. Built to delineate standing dead trees
from a hand-digitised point layer — the operator supplies the objects, the
algorithm supplies their boundaries.

```r
grow_seeds_raster("ortho_lab.tif", "dead_trees.shp",
                  output = "labels.tif", polygons = "crowns.gpkg",
                  max_cost = 15, band_weights = c(0.5, 2.5, 1),
                  max_radius = 20, fill_holes = TRUE)
```

`max_cost` is a tolerance in the band units — a ΔE tolerance when the input is
CIELAB, so feed Lab, not RGB. `band_weights` reshapes the feature space
(weighting `a*` up separates dead brown from living green); `max_radius` bounds
the reach; `fill_holes` closes the pockets a cut leaves inside a crown. Label
`i` is the region grown from the i-th point, so it joins back to that point's
attributes. `docs/grow_seeds_guide.md` is the operator's guide, with a worked
recipe for dead trees.

It runs on the same frozen IFT kernel as `sicle`, which is what keeps it
bit-identical to the Python twin.

## Agreement with plGeoAdaptels

Bit-identical, and checked rather than asserted:

```sh
pip install plgeoadaptels
python3 tools/generate_plgeoadaptels_reference.py
Rscript tools/cross_validate_against_plgeoadaptels.R
```

Thirty cases across all three algorithms — every adaptels metric and both
connectivities, a non-default Minkowski exponent, the normalise path, a mask
with an interior hole, single- and multi-band input, a constant raster,
`enforce_connectivity` at three `min_size` values, SICLE across its parameters
and a saliency map, and `grow_seeds` across every option including the
`fill_holes` cleanup. Zero differing pixels in all of them. The real 400x400
three-band scene in the plGeoAdaptels repository — and the CIELAB dead-tree
ortho with 36 hand-placed points — were checked the same way locally; they are
not shipped here.

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
Raster Data*. R package version 0.3.0.
https://github.com/igorpawelec/rgeoadaptels

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

## Licence

GPL-3.
