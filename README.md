# rgeoadaptels

<img src="https://raw.githubusercontent.com/igorpawelec/rgeoadaptels/main/www/rgeoadaptels.png" alt="rgeoadaptels logo" align="right" width="200"/>

[![R-CMD-check](https://github.com/igorpawelec/rgeoadaptels/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igorpawelec/rgeoadaptels/actions/workflows/R-CMD-check.yaml)
[![Release](https://img.shields.io/github/v/release/igorpawelec/rgeoadaptels)](https://github.com/igorpawelec/rgeoadaptels/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%203.6-blue.svg)](https://www.r-project.org)

**Scale-adaptive superpixels (adaptels) and SICLE superpixels for geospatial raster data — and `grow_seeds`, the same kernel run the other way round.**

A raster goes in; a segmentation comes out, as a label matrix or a raster. R and plain C — no Rcpp, no required packages, and `terra` only if you want to read and write files.

> **Python users:** the same algorithms are in [pygeoadaptels](https://github.com/igorpawelec/pygeoadaptels). The two are separate repositories because their tooling and idioms do not mix, but they are **bit-identical** — see [Agreement](#agreement-with-pygeoadaptels).

## The problem it solves

Object-based analysis of an orthophoto starts by cutting the image into segments, and the usual superpixel methods ask you how many. A fixed count imposes a grid on a scene that has no grid: the same segment size that resolves a small crown is wasted on a stretch of shadow, and the count that suits one plot is wrong on the next.

**Adaptels** let the scene decide. A region grows from a seed until its internal colour distance passes a threshold, and the pixels beyond it become seeds in turn. Where the image is textured the regions stay small, where it is homogeneous they grow large, and the count follows from the data. One parameter, and it is a distance in the units of your bands.

**`grow_seeds`** is the inverse. When you already know where the objects are — a point layer of standing dead trees digitised by an operator — you do not want a partition of the whole image, you want the boundary of each object. Every point grows into the region that looks like the pixel it sits on, within a tolerance, and everything unseeded stays unassigned. The operator supplies the objects, the algorithm supplies their extent.

Both run on one frozen kernel — an image foresting transform with a min-heap — which is also what **SICLE** uses when you do want a fixed number of superpixels.

<img src="https://raw.githubusercontent.com/igorpawelec/rgeoadaptels/main/www/pipeline.png" alt="A 0.25 m orthophoto of a circular spruce plot with 36 dead-tree points; the same plot partitioned into 2,792 adaptels at threshold 60; and the crowns grown from the 36 points by grow_seeds" width="100%"/>

*A 100 × 100 m circular sample plot in a spruce stand at 0.25 m (`test_data/SNP_21_2020_1.tif`), with 36 standing dead trees digitised as points. Middle: the plot partitioned into 2,792 adaptels at the default threshold of 60 — median 1.8 m², largest 50 m². Right: the crowns `grow_seeds` grows from the 36 points on the CIELAB version of the same scene, with the dead-tree recipe (`max_cost` 15, `band_weights` 0.5/2.5/1, `max_radius` 20, `fill_holes`). The orthophoto is contrast-stretched for display only. Made by `www/figures.R`, from this package — and every number in these captions is the same as in pygeoadaptels' README, to the adaptel and the pixel.*

## The one parameter that matters

`threshold` is the colour distance a region may accumulate before it stops growing, in the units of the input bands. It is not a size and not a count: it decides the *scale* at which the scene is cut, and the scene then decides how many pieces that takes.

<img src="https://raw.githubusercontent.com/igorpawelec/rgeoadaptels/main/www/threshold.png" alt="A 30 by 30 m window of the plot at threshold 30, 60 and 120: 825, 403 and 178 adaptels, small on the textured crowns and large in the shadow" width="100%"/>

*A 30 × 30 m window of the plot at threshold 30, 60 and 120 — 825, 403 and 178 adaptels in the window, 6,138, 2,792 and 1,159 on the whole plot. At every threshold the adaptels are small on the textured crowns and large in the homogeneous shadow between them; the threshold shifts that whole distribution rather than fixing a size.*

The threshold is **per metric, not universal**, and this is the single easiest way to get a nonsense result, so the package refuses rather than obliges:

| distance | range | sensible threshold |
|---|---|---|
| `minkowski` | grows with the data range | 60 on 0–255 imagery |
| `cosine` | [0, 1] by construction | around 0.03 |
| `angular` | [0, 1] by construction | around 0.03 |

Passing 60 to `cosine` would merge the whole raster into one adaptel. It raises instead, and says what to try. The same applies to `normalize = TRUE`: normalising caps the largest possible Minkowski distance at `n_bands^(1/p)` — about 1.73 for three bands — so the raw-data default cannot be reached and is rejected. `cosine` and `angular` compare the *direction* of the spectral vector rather than its length, so the same material lit differently — a crown in sun versus in shade — lands in one adaptel; `minkowski` will split it.

## The other way round: `grow_seeds`

`grow_seeds` takes the points and one tolerance, `max_cost`: how far from the seed pixel's colour a crown may reach, in band units. Feed it CIELAB — convert with [rgeopalette](https://github.com/igorpawelec/rgeopalette) first — and `max_cost` becomes a ΔE, a colour difference with a meaning.

<img src="https://raw.githubusercontent.com/igorpawelec/rgeoadaptels/main/www/max_cost.png" alt="The 36 dead-tree points grown at max_cost 8, 15 and 25: crowns covering 1.9, 4.9 and 5.5 percent of the plot; at 25 one crown floods to a disc that max_radius stops" width="100%"/>

*The same 36 points grown at `max_cost` 8, 15 and 25: crowns covering 1.9 %, 4.9 % and 5.5 % of the plot, median 2.5, 6.4 and 6.5 m². At 8 the crowns stop short of their own edges; at 15, the recipe, they fill the bleached crowns and stop at the living neighbours; at 25 one seed on a dark stem leaks into the surrounding canopy until `max_radius` (20 px, 5 m) stops it — the disc in the lower left is that cap doing its job. Calibrate `max_cost` by sweeping it and looking, exactly like this.*

`band_weights` reshapes the feature space — weighting `a*` up separates dead brown from living green, which is the whole dead-vs-living problem — `max_radius` bounds the reach, and `fill_holes` closes the pockets a cut leaves inside a crown. Label `i` is the region grown from the i-th point, so it joins back to that point's attributes. `docs/grow_seeds_guide.md` is the operator's guide, with the worked recipe for dead trees.

## When to use it, and when not

Use adaptels when the segments are the units of a later analysis — zonal statistics, a classifier over segment features, a manual interpretation — and you want them to follow the scene rather than a grid. Use `grow_seeds` when the objects are already located and you need their boundaries. Use SICLE when a downstream method needs a fixed number of superpixels, or a saliency map (a normalised height model, say) should pull the boundaries towards objects.

Adaptels are not object detection: a segment is a homogeneous patch, not a tree, and a crown may be several of them. For crowns from a canopy height model, use [rcacumen](https://github.com/igorpawelec/rcacumen).

### The package family

rgeoadaptels is one step of a longer chain; the other steps are separate packages, each with an R and a Python twin.

| Step | R | Python |
|---|---|---|
| Colour-space conversion of orthophotos | [rgeopalette](https://github.com/igorpawelec/rgeopalette) | [pygeopalette](https://github.com/igorpawelec/pygeopalette) |
| Adaptive superpixels and seeded growing on orthophotos | **rgeoadaptels** | [pygeoadaptels](https://github.com/igorpawelec/pygeoadaptels) |
| Crowns from a canopy height model | [rcacumen](https://github.com/igorpawelec/rcacumen) | [pycacumen](https://github.com/igorpawelec/pycacumen) |
| Standing dead trees on orthophotos | — | [pygeosnag](https://github.com/igorpawelec/pygeosnag) |
| The same, inside QGIS | | [qgis-geoadaptels-geopalette](https://github.com/igorpawelec/qgis-geoadaptels-geopalette), [qgis-geosnag](https://github.com/igorpawelec/qgis-geosnag) |
| Polish national geodata (GUGiK, BDL) | [rgeopl](https://github.com/igorpawelec/rgeopl) | — |

## Install

```r
# install.packages("remotes")
remotes::install_github("igorpawelec/rgeoadaptels")
```

A C compiler is needed, which R already requires on Linux and macOS; on Windows install [Rtools](https://cran.r-project.org/bin/windows/Rtools/). `terra` is optional and only used to read and write rasters.

## Use

```r
library(rgeoadaptels)

# data is (bands, rows, cols), or a plain matrix for one band
seg <- adaptels(data, threshold = 60)
seg$n_adaptels
seg$labels          # integer matrix, 0-based ids, -9999 for nodata

adaptels_raster("scene.tif", "adaptels.tif", threshold = 60)   # with a raster, via terra
```

```r
grow_seeds_raster("ortho_lab.tif", "dead_trees.shp",
                  output = "labels.tif", polygons = "crowns.gpkg",
                  max_cost = 15, band_weights = c(0.5, 2.5, 1),
                  max_radius = 20, fill_holes = TRUE)
```

```r
seg <- sicle(data, n_segments = 200)
seg$labels            # 0-based ids, each a single 8-connected region
seg$n_superpixels
sicle(data, seeds = cbind(rows, cols), n_segments = 200)   # seeds given instead of sampled; 1-based (row, col)
```

## Reference

### Adaptels can arrive in more than one piece

Adaptels compete: a later one takes a pixel from an earlier one whenever it arrives with a smaller accumulated distance. That competition is what gives the method its boundary adherence, and it can also cut an earlier adaptel in two. On a 400 × 400 three-band scene at the default threshold, about 10 per cent come out split. Harmless if the labels are a lookup; not harmless for zonal statistics, which would average two spatially separate patches into one "object":

```r
seg   <- adaptels(data, threshold = 60)
split <- enforce_connectivity(seg$labels)     # every region now contiguous
```

Not applied automatically, because it changes the adaptel count.

### `grow_seeds`

`grow_seeds(data, seeds, mask = NULL, max_cost = NULL, band_weights = NULL, compactness = 0, seed_window = 1L, max_radius = NULL, fill_holes = FALSE)`. Seeds are a `(n, 2)` matrix of 1-based `(row, col)`; `grow_seeds_raster` takes a point layer in any CRS and does the conversion. Labels are 0-based, `-1` unassigned, and label `i` is the region grown from `seeds[i, ]`. `max_cost = NULL` keeps every reachable pixel — a partition, as the kernel does for SICLE.

### SICLE

`sicle` starts from far more seeds than you want (`n_oversampling`, 3000), grows an optimum-path forest, scores every seed and discards the least relevant, and repeats until `n_segments` remain. Seeds can be given instead of sampled, which is what lets the two twins be compared on the algorithm rather than the sampler — NumPy's `Generator.choice` cannot be reproduced outside NumPy. `n_iterations` does less than it looks: the paper's preservation curve makes 3 bit-identical to 2. See `?sicle` for why, and why 2 is a speed setting rather than a quality one.

### Agreement with pygeoadaptels

Bit-identical, and checked rather than asserted:

```sh
pip install pygeoadaptels
python3 tools/generate_pygeoadaptels_reference.py
Rscript tools/cross_validate_against_pygeoadaptels.R
```

Thirty cases across all three algorithms — every adaptels metric and both connectivities, a non-default Minkowski exponent, the normalise path, a mask with an interior hole, single- and multi-band input, a constant raster, `enforce_connectivity` at three `min_size` values, SICLE across its parameters and a saliency map, and `grow_seeds` across every option including the `fill_holes` cleanup. Zero differing pixels in all of them. The figures above are a fourth check in plain sight: `www/figures.R` and pygeoadaptels' `www/figures.py` run the same scene through the same parameters and print the same counts, sizes and coverages.

<details>
<summary><b>Why an equality is possible here, and not in the other twin pairs</b></summary>

The check is an **equality**, not a tolerance, and that is worth a note because the other twin pairs in this family could not manage it. rcacumen reimplements `scikit-image`'s watershed and differs on 0.25 % of pixels where plateau ties fall differently. rgeopalette computes in double where pygeopalette stores single precision. Neither problem exists here: this is a port of the same kernel — same heap, same neighbour order, same arithmetic — so there is no second implementation to disagree with, and anything short of identical would be a bug.

Three details make that possible, and all three are easy to get wrong:

- **R stores matrices by column, Python by row.** The C is handed `t(matrix)`, whose column-major buffer is byte-identical to a row-major NumPy array. Skip it and the neighbour order transposes, which changes which adaptel claims a contested pixel.
- **Bands are indexed `layers[l * size + i]`**, matching NumPy's `(n_layers, size)` C-order rather than an R matrix's layout.
- **The heap is 1-based**, as in the Python, so the sift arithmetic is the same expression rather than a translated one.

</details>

<details>
<summary><b>The figures</b></summary>

`Rscript www/figures.R` from the repository root remakes the three figures from the files in `test_data/` — the orthophoto, its CIELAB version made with pygeopalette, and the 36 dead-tree points — using the package in the checkout (`pkgload`) rather than an installed copy, and prints the numbers the captions quote.

</details>

## Citation

Pawelec, I. (2026). *rgeoadaptels: Scale-Adaptive Superpixels for Geospatial Raster Data*. R package version 0.4.0. https://github.com/igorpawelec/rgeoadaptels

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff), including the papers the algorithms come from: Achanta et al. (2018) for adaptels and Belém et al. (2023) for SICLE. The original C implementation of adaptels, `plGeoAdaptels`, was written by Paweł Netzel at the University of Agriculture in Kraków.

## Licence

GPL-3.
