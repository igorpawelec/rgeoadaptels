# Prove rgeoadaptels and plGeoAdaptels produce the same segmentation.
#
# This is an equality check, not a tolerance, and deliberately so. The twin
# packages in this family had to settle for tolerances: rHRG reimplements
# scikit-image's watershed and differs on 0.25% of pixels where plateau ties
# fall differently, and GeoPaletteR computes in double where GeoPalette
# stores single precision. Neither applies here. This is a port of the same
# kernel — same heap, same neighbour order, same arithmetic — with no second
# implementation to disagree with, so anything short of bit-identical is a
# bug rather than a difference of opinion.
#
#   python3 tools/generate_plgeoadaptels_reference.py
#   Rscript tools/cross_validate_against_plgeoadaptels.R
#
# Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.

library(rgeoadaptels)

dir <- file.path("tools", "reference")
if (!dir.exists(dir))
  stop("run tools/generate_plgeoadaptels_reference.py first", call. = FALSE)

rd <- function(f) as.matrix(utils::read.csv(file.path(dir, f), header = FALSE))

scene <- function(name, nb) {
  b <- lapply(seq_len(nb) - 1L,
              function(l) rd(sprintf("%s_band%d.csv", name, l)))
  a <- array(0, dim = c(nb, nrow(b[[1]]), ncol(b[[1]])))
  for (l in seq_len(nb)) a[l, , ] <- b[[l]]
  a
}
SCENES <- list(multi = scene("multi", 3L), single = scene("single", 1L),
               flat = scene("flat", 2L))
mask <- rd("mask.csv")

CASES <- list(
  mink_rook  = quote(adaptels(d, threshold = 60)),
  mink_queen = quote(adaptels(d, threshold = 60, queen_topology = TRUE)),
  mink_p3    = quote(adaptels(d, threshold = 60, minkowski_p = 3)),
  cosine     = quote(adaptels(d, threshold = 0.03, distance = "cosine")),
  angular    = quote(adaptels(d, threshold = 0.03, distance = "angular")),
  normalize  = quote(adaptels(d, threshold = 0.4, normalize = TRUE)),
  masked     = quote(adaptels(d, mask = mask, threshold = 60)),
  single     = quote(adaptels(d, threshold = 60)),
  flat       = quote(adaptels(d, threshold = 60)),
  tight      = quote(adaptels(d, threshold = 8))
)
SCENE_OF <- c(mink_rook = "multi", mink_queen = "multi", mink_p3 = "multi",
              cosine = "multi", angular = "multi", normalize = "multi",
              masked = "multi", single = "single", flat = "flat",
              tight = "multi")

cases <- utils::read.csv(file.path(dir, "cases.csv"), header = FALSE,
                         col.names = c("tag", "scene", "n"),
                         colClasses = c("character", "character", "integer"))

bad <- 0L
cat(sprintf("%-11s %8s %8s %10s %10s\n", "case", "R", "python", "identical",
            "diff px"))
cat(strrep("-", 52), "\n")

for (tag in names(CASES)) {
  d <- SCENES[[SCENE_OF[[tag]]]]
  got <- eval(CASES[[tag]])
  ref <- rd(sprintf("out_%s.csv", tag))
  n_py <- cases$n[cases$tag == tag]
  same <- identical(as.vector(got$labels), as.vector(ref)) &&
    identical(got$n_adaptels, n_py)
  ndiff <- sum(got$labels != ref)
  if (!same) bad <- bad + 1L
  cat(sprintf("%-11s %8d %8d %10s %10d\n", tag, got$n_adaptels, n_py,
              same, ndiff))
}

# ── SICLE ────────────────────────────────────────────────────────────
#
# Seeds come from the reference so that the sampler is out of the
# comparison: NumPy's Generator.choice cannot be reproduced outside NumPy,
# and Belem et al. treat the sampling as a free choice rather than part of
# the algorithm. The +1 is the 0-based to 1-based crossing.
sicle_seeds <- rd("sicle_seeds.csv") + 1L
sicle_sal <- rd("sicle_saliency.csv")
d <- SCENES[["multi"]]

SICLE_CASES <- list(
  s_it2  = quote(sicle(d, seeds = sicle_seeds, n_segments = 60,
                       n_iterations = 2, quiet = TRUE)),
  s_it5  = quote(sicle(d, seeds = sicle_seeds, n_segments = 60,
                       n_iterations = 5, quiet = TRUE)),
  s_it10 = quote(sicle(d, seeds = sicle_seeds, n_segments = 60,
                       n_iterations = 10, quiet = TRUE)),
  s_many = quote(sicle(d, seeds = sicle_seeds, n_segments = 200,
                       n_iterations = 2, quiet = TRUE)),
  s_few  = quote(sicle(d, seeds = sicle_seeds, n_segments = 5,
                       n_iterations = 2, quiet = TRUE)),
  s_sal  = quote(sicle(d, seeds = sicle_seeds, n_segments = 60,
                       n_iterations = 2, saliency = sicle_sal, quiet = TRUE)),
  s_mask = quote(sicle(d, seeds = rd("sicle_seeds_masked.csv") + 1L,
                       mask = mask, n_segments = 60, n_iterations = 2,
                       quiet = TRUE))
)

for (tag in names(SICLE_CASES)) {
  got <- eval(SICLE_CASES[[tag]])
  ref <- rd(sprintf("out_%s.csv", tag))
  n_py <- cases$n[cases$tag == tag]
  same <- identical(as.vector(got$labels), as.vector(ref)) &&
    identical(got$n_superpixels, n_py)
  ndiff <- sum(got$labels != ref)
  if (!same) bad <- bad + 1L
  cat(sprintf("%-11s %8d %8d %10s %10d
", tag, got$n_superpixels,
              n_py, same, ndiff))
}

base <- rd("conn_input.csv")
for (ms in c(0, 5, 20)) {
  got <- enforce_connectivity(base, min_size = ms)
  ref <- rd(sprintf("conn_%d.csv", ms))
  n_py <- cases$n[cases$tag == sprintf("conn%d", ms)]
  same <- identical(as.vector(got$labels), as.vector(ref)) &&
    identical(got$n_adaptels, n_py)
  ndiff <- sum(got$labels != ref)
  if (!same) bad <- bad + 1L
  cat(sprintf("%-11s %8d %8d %10s %10d\n", sprintf("conn ms=%d", ms),
              got$n_adaptels, n_py, same, ndiff))
}

cat(strrep("-", 52), "\n")
if (bad > 0L)
  stop(bad, " case(s) differ from plGeoAdaptels", call. = FALSE)
cat("rgeoadaptels == plGeoAdaptels, bit for bit, on all ",
    length(CASES) + length(SICLE_CASES) + 3L, " cases.\n", sep = "")
