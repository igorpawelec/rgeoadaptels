# Plain-R tests: any error fails R CMD check. No testthat dependency, so the
# package stays free of test-only requirements.
#
# These do not compare against plGeoAdaptels — that needs Python and lives in
# tools/. They pin the behaviour this package promises on its own.

library(rgeoadaptels)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

scene <- function(bands = 3L, rows = 30L, cols = 41L, seed = 1L) {
  set.seed(seed)
  array(as.numeric(sample(0:255, bands * rows * cols, TRUE)),
        dim = c(bands, rows, cols))
}

# ── shape and structure ──────────────────────────────────────────────
cat("structure\n")

d <- scene()
o <- adaptels(d, threshold = 60)
ok(identical(dim(o$labels), c(30L, 41L)), "labels keep the raster's shape")
ok(is.integer(o$labels), "labels are integer")
ok(o$n_adaptels > 1L, "a noise scene gives a real segmentation")
ok(identical(sort(unique(as.vector(o$labels))), 0:(o$n_adaptels - 1L)),
   "ids are consecutive from 0 with no gaps")

# Non-square on purpose: a transposed buffer would still fit a square scene.
ok(nrow(o$labels) != ncol(o$labels), "the test scene is not square")

o2 <- adaptels(matrix(as.numeric(sample(0:255, 26 * 19, TRUE)), 26, 19),
               threshold = 60)
ok(identical(dim(o2$labels), c(26L, 19L)), "a 2-D matrix is accepted")

# ── known behaviour ──────────────────────────────────────────────────
cat("behaviour\n")

flat <- array(42, dim = c(2L, 20L, 26L))
ok(adaptels(flat, threshold = 60)$n_adaptels == 1L,
   "a constant raster is one adaptel")

lo <- adaptels(d, threshold = 10)$n_adaptels
hi <- adaptels(d, threshold = 300)$n_adaptels
ok(lo > hi, "a lower threshold gives more adaptels")

n4 <- adaptels(d, threshold = 60)$n_adaptels
n8 <- adaptels(d, threshold = 60, queen_topology = TRUE)$n_adaptels
ok(n8 <= n4, "8-connectivity gives no more adaptels than 4")

# Each metric on its own scale. `n >= 1` would pass even when a metric
# collapses the raster into a single adaptel, which is how an inverted
# cosine went unnoticed in plGeoAdaptels for two releases.
for (spec in list(c("minkowski", "60"), c("cosine", "0.03"),
                  c("angular", "0.03"))) {
  n <- adaptels(d, threshold = as.numeric(spec[2]), distance = spec[1])$n_adaptels
  ok(n > 1L, paste0("'", spec[1], "' produces a real segmentation, not one adaptel"))
}

a <- adaptels(d, threshold = 60)
b <- adaptels(d, threshold = 60)
ok(identical(a$labels, b$labels), "repeated runs are identical")

# ── nodata ───────────────────────────────────────────────────────────
cat("nodata\n")

m <- matrix(0L, 30, 41)
m[1, ] <- 1L; m[, 41] <- 1L; m[10:14, 15:20] <- 1L
om <- adaptels(d, mask = m, threshold = 60)
ok(all(om$labels[m == 1] == -9999), "masked pixels come back as -9999")
ok(all(om$labels[m == 0] >= 0), "unmasked pixels are all labelled")
ok(identical(sort(unique(as.vector(om$labels[m == 0]))),
             0:(om$n_adaptels - 1L)),
   "ids stay consecutive when a mask is present")

dna <- d
dna[1, 5, 5] <- NA
ok(adaptels(dna, threshold = 60)$labels[5, 5] == -9999,
   "NA in a band is treated as nodata when no mask is given")

# ── enforce_connectivity ─────────────────────────────────────────────
cat("connectivity\n")

# Grown with queen_topology so that adaptels exist which are 8-connected but
# not 4-connected; enforce_connectivity tests 4, so the splitting path runs.
# The default rook growth splits nothing on noise and would prove nothing.
q <- adaptels(d, threshold = 200, queen_topology = TRUE)
s0 <- enforce_connectivity(q$labels)
ok(s0$n_adaptels > q$n_adaptels, "8-connected growth splits under a 4-connected check")
ok(all((q$labels < 0) == (s0$labels < 0)), "nodata stays nodata through the split")

s5 <- enforce_connectivity(q$labels, min_size = 5)
ok(s5$n_adaptels < s0$n_adaptels, "min_size absorbs slivers")

r <- adaptels(d, threshold = 60)
same <- enforce_connectivity(r$labels)
ok(same$n_adaptels >= r$n_adaptels, "splitting never reduces the count at min_size 0")

# A single-component adaptel is kept whole however small, because the
# min_size rule only applies once an adaptel has come apart.
tiny <- matrix(-9999L, 5, 5)
tiny[3, 3] <- 0L
t1 <- enforce_connectivity(tiny, min_size = 100)
ok(t1$n_adaptels == 1L && t1$labels[3, 3] == 0L,
   "a lone one-pixel adaptel survives a large min_size")

hand <- matrix(-9999L, 6, 6)
hand[2, 2] <- 0L; hand[5, 5] <- 0L          # same id, nowhere near
hand[2, 5] <- 1L
h <- enforce_connectivity(hand)
ok(h$n_adaptels == 3L, "a hand-built split becomes two adaptels plus the other")
ok(h$labels[2, 2] != h$labels[5, 5], "the two pieces get different ids")

# A raster with nothing in it returns -1, not -9999, and that is deliberate
# rather than an oversight here: plGeoAdaptels takes an early return in this
# case which skips the line that restores the original nodata marker, so a
# raster holding one valid pixel comes back with -9999 while a raster holding
# none comes back with -1. The inconsistency is upstream; replicating it is
# what keeps the cross-check in tools/ an equality rather than a tolerance.
allnod <- matrix(-9999L, 4, 4)
an <- enforce_connectivity(allnod)
ok(an$n_adaptels == 0L && all(an$labels == -1L),
   "an empty raster returns -1, matching plGeoAdaptels' early return")

one <- matrix(-9999L, 4, 4); one[2, 2] <- 0L
on <- enforce_connectivity(one)
ok(on$n_adaptels == 1L && on$labels[1, 1] == -9999L,
   "one valid pixel is enough for nodata to keep its -9999")

# ── errors ───────────────────────────────────────────────────────────
cat("validation\n")

err <- function(expr, what) {
  ok(inherits(try(expr, silent = TRUE), "try-error"), what)
}
err(adaptels(d, threshold = -1), "a negative threshold is rejected")
err(adaptels(d, threshold = 0), "a zero threshold is rejected")
err(adaptels(d, threshold = 60, distance = "cosin"),
    "a misspelt distance is rejected rather than falling back")
err(adaptels(d, threshold = 60, distance = "cosine"),
    "a minkowski-scaled threshold is rejected for cosine")
err(adaptels(d, threshold = 60, normalize = TRUE),
    "the raw-data default is rejected once bands are normalised")
err(adaptels("not a raster", threshold = 60), "non-numeric data is rejected")
err(adaptels(array(0, dim = c(2, 2, 2, 2)), threshold = 60),
    "a 4-D array is rejected")
err(adaptels(d, mask = matrix(0, 3, 3), threshold = 60),
    "a mask of the wrong size is rejected")
err(enforce_connectivity(array(0, dim = c(2, 2, 2))),
    "enforce_connectivity rejects a 3-D input")
err(enforce_connectivity(matrix(0L, 4, 4), min_size = -1),
    "a negative min_size is rejected")

# The normalised scale is usable, not merely guarded.
ok(adaptels(d, threshold = 0.4, normalize = TRUE)$n_adaptels > 1L,
   "a threshold on the normalised scale works")

cat("\nall adaptel tests passed\n")
