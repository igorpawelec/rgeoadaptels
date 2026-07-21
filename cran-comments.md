# cran-comments

## Test environments

* Windows 11, R 4.5.2, Rtools45 (local), `R CMD check --as-cran`, with
  TinyTeX installed so that `checking PDF version of manual` runs rather
  than being skipped
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)

## R CMD check results

0 errors, 0 warnings, 1 note.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Igor Pawelec <igor.pawelec@urk.edu.pl>'
New submission
```

Expected for a first submission.

The remaining local notes are environmental rather than properties of the
package: `checking top-level files` reports that README.md and NEWS.md
cannot be checked without pandoc, and `checking for future file
timestamps` reports that it could not reach the network time service.
Neither appears on the CI runners, which install pandoc.

## Notes for the reviewer

The package compiles two C files and uses no Rcpp: the region grower,
which is a sequential priority flood, and a connected-component pass.
Everything else is plain R. There are no imports beyond base R; `terra`
is a Suggests used only by `read_bands()` and `adaptels_raster()`, and
the segmentation itself runs on ordinary numeric arrays.

`tools/` is excluded from the build via `.Rbuildignore`. It holds a
cross-check against the plGeoAdaptels Python package, which needs Python
and so cannot run on a CRAN machine.
