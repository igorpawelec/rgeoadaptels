# Contributing to rgeoadaptels

rgeoadaptels segments rasters into adaptels (Achanta 2018) and SICLE
superpixels (Belém 2023). It is the R twin of
[plGeoAdaptels](https://github.com/igorpawelec/plGeoAdaptels); the two are
developed together and are meant to agree.

Bug reports, ideas and pull requests are all welcome.

## Reporting a bug

Please include:

- what you ran, ideally a snippet that reproduces it,
- what you expected and what happened,
- the raster's dimensions, type and value range — most surprises trace back
  to nodata, units or an unexpected range,
- `sessionInfo()`.

A synthetic input that reproduces the problem helps more than a real one.

## Setup

```r
install.packages(c("roxygen2", "devtools"))
devtools::install_deps(dependencies = TRUE)
devtools::load_all()
```

The algorithms are plain C in `src/`, called through `.Call` — no Rcpp.
`src/heap.h` holds the 1-based growable heap shared by the adaptel and SICLE
code. Symbols are registered in `src/init.c` with `R_forceSymbols(dll, TRUE)`,
so a new `.Call` entry point has to be added there or it will not be found.

Matrices cross the boundary transposed. R is column-major and NumPy is
row-major, so `t(m)` is what makes the two buffers line up; getting this
wrong produces a result that looks plausible and is transposed.

## Running the checks

```bash
R CMD build .
R CMD check --as-cran rgeoadaptels_X.Y.Z.tar.gz
```

Tests are plain R scripts in `tests/`, run by `R CMD check`; a failure there
is an ERROR, not a note. Run the full check rather than the tests alone — the
things that have actually broken this family were documentation and metadata
problems the tests cannot see.

## Documentation

`man/` is generated. Edit the roxygen comments above the function in `R/`,
then run `roxygen2::roxygenise()`; never edit a `.Rd` by hand.

This is not a style preference. Every file in rHRG's `man/` once carried a
hand-written header instead of roxygen2's marker, so `roxygenise()` silently
skipped all of them and the directory drifted from the comments it claimed to
come from. The only part `R CMD check` notices is argument defaults, via
`codoc`; prose can be wrong indefinitely.

## Agreement with plGeoAdaptels

```bash
python3 tools/generate_plgeoadaptels_reference.py   # needs plgeoadaptels installed
Rscript tools/cross_validate_against_plgeoadaptels.R
```

This is the one twin pair that is **bit-identical**, and it is worth keeping
that way: both sides take explicit seeds, so there is no RNG for the two to
diverge on. Any difference at all is a defect in one of them.

**Run it if you change the growing, the connectivity or the SICLE
reduction.**

Never compare two segmentations by label equality. Region numbering is an
artefact of iteration order, so a renumbering reads as a 98 % difference when
nothing moved. Match regions by best overlap first; this has produced a false
alarm in this family three separate times.

## The pinned twin

The agreement job installs plGeoAdaptels from a **tag**, not the tip of its main
branch — see `.github/workflows/R-CMD-check.yaml`. Unpinned, a commit in
plGeoAdaptels could turn this repository's CI red without anything changing here,
and the failure would read as a defect in this package.

The cost is that the pin goes stale. When plGeoAdaptels releases and the segmentation
changes, bump the tag in the workflow and re-run the agreement job. A stale
pin proves agreement against an old plGeoAdaptels, which is quieter than a red
build and therefore worse.

## Pull requests

- Add a test that fails before your change and passes after. Confirm it fails
  — a test written after the fix and never seen red guards nothing. Several
  in this family were caught doing exactly that, by sabotaging the code the
  test was supposed to protect.
- If you change what the package produces, say so in the pull request and in
  `NEWS.md`. Silent changes to output are the hardest kind to debug for
  anyone with a pipeline in flight.
- Numbers in documentation should come from a measurement, not an estimate.

## Releasing

The checklist exists because of a specific failure. `max_iters` changed
default in 0.3.0, and that one change broke CI in both HRG twins at once —
pyHRG with `int | None`, a runtime `TypeError` before Python 3.10, and rHRG
with a `man/` page still documenting `200L`, which `R CMD check` reports as a
codoc WARNING and the workflow treats as a failure. Neither was noticed.
**pyHRG then tagged 0.3.0, 0.4.0 and 0.5.0 with the workflow red**; rHRG
shipped two the same way and rgeoadaptels two more.

Local checks passed in every one of those cases. They were run on one R
version, on one operating system, by someone who already knew what the change
was meant to do. The matrix is the part that disagrees.

1. Update `NEWS.md`. If the output changes, say so in those words.
2. **`roxygen2::roxygenise()`.** Skipping this is what broke rHRG: `man/`
   kept a default the code no longer had, and `?grow_crowns` told readers the
   cap was still there.
3. Bump `Version:` in `DESCRIPTION`, and the version and `date-released` in
   `CITATION.cff`. Search for the *old* number and read the hits —
   `grep -rn "0.1.0" --exclude-dir=.git` — rather than editing the places you
   remember.
4. `R CMD build .` then `R CMD check --as-cran <tarball>`. It must end
   `0 errors | 0 warnings`; the workflow fails on a warning, so a warning
   here is a red build there. Two NOTEs are expected offline: missing pandoc,
   and new-submission.
5. Commit and push. **Do not tag yet.**
6. **Wait for Actions on the pushed commit and confirm every matrix job is
   green.** Not the previous run, not the branch generally — that commit.
   This is the step that was missing. Either open the Actions tab, or:

   ```bash
   curl -s "https://api.github.com/repos/igorpawelec/rgeoadaptels/actions/runs?per_page=1" | python -c "import json,sys; r=json.load(sys.stdin)['workflow_runs'][0]; print(r['head_sha'][:7], r['status'], r['conclusion'])"
   ```

   `gh run list` is nicer if the GitHub CLI is installed; it is not
   everywhere, and the curl form needs nothing but a public repo.
7. Only then tag and push the tag:
   `git tag -a vX.Y.Z -m "..." && git push --tags`

The order matters. A tag is what people install and what a DOI points at, so
it should never be the thing that discovers a broken build. If Actions is
red, fix it and release the fix as its own version — the broken tag stays in
history either way.

## Licence

rgeoadaptels is GPLv3. Contributions are accepted under the same licence.
