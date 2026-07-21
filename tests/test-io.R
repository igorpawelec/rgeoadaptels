# The terra bridge — the one place in the package where the data changes
# shape, and therefore the one place a silent transpose can hide.
#
# terra is a Suggests, so this file is a no-op when it is absent.

library(rgeoadaptels)

if (!requireNamespace("terra", quietly = TRUE)) {
  cat("terra not installed; skipping the raster bridge tests\n")
} else {

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

# Not square, and each cell encodes its own position. A 4x4 of random values
# cannot tell a correct write from a transposed one: the dimensions still
# match and nothing looks wrong.
NR <- 9L; NC <- 14L
pos <- outer(seq_len(NR), seq_len(NC), function(r, c) r * 10 + c)

mk <- function(nl = 3L, vals = NULL) {
  r <- terra::rast(nrows = NR, ncols = NC, nlyrs = nl,
                   xmin = 0, xmax = NC, ymin = 0, ymax = NR)
  if (is.null(vals))
    vals <- do.call(cbind, lapply(seq_len(nl),
                                  function(l) as.vector(t(pos)) + (l - 1) * 100))
  terra::values(r) <- vals
  r
}

cat("read_bands\n")

b <- read_bands(mk())
ok(identical(dim(b$data), c(3L, NR, NC)), "bands come back as (bands, rows, cols)")
ok(identical(b$data[1, , ], pos), "band 1 values land in the right cells")
ok(identical(b$data[2, , ], pos + 100), "band 2 is the second layer")
ok(all(b$mask == 0L), "a raster without NA has an empty mask")
ok(inherits(b$template, "SpatRaster"), "the template comes back")

r_na <- mk()
v <- terra::values(r_na); v[5, 1] <- NA
terra::values(r_na) <- v
ok(sum(read_bands(r_na)$mask) == 1L, "NA in a band becomes one masked pixel")

cat("adaptels_raster\n")

set.seed(4)
rr <- mk(3L, matrix(as.numeric(sample(0:255, NR * NC * 3, TRUE)), ncol = 3))
seg <- adaptels_raster(rr, output = NULL, threshold = 60)
ok(seg$n_adaptels > 1L, "a noise raster segments into more than one adaptel")
ok(identical(dim(seg$labels), c(NR, NC)), "labels keep the raster's shape")

b2 <- read_bands(rr)
direct <- adaptels(b2$data, mask = b2$mask, threshold = 60)
ok(identical(seg$labels, direct$labels),
   "adaptels_raster agrees with adaptels on the same bands")

# The round trip. Write, read back with a fresh terra call, require the
# values to be where they were. A missing t() looks plausible on a square
# raster; 9x14 is where it stops looking plausible.
f <- tempfile(fileext = ".tif")
on.exit(unlink(f), add = TRUE)
invisible(adaptels_raster(rr, output = f, threshold = 60))
ok(file.exists(f), "an output file is written")

back <- terra::rast(f)
got <- terra::as.matrix(back, wide = TRUE)
# The shape check alone proves nothing: terra keeps the raster geometry
# whatever order the values arrive in, so a transposed write still returns
# 9x14. Only the value comparison catches it.
ok(identical(dim(got), c(NR, NC)), "the output keeps its 9x14 geometry")
ok(all(got == seg$labels), "values round-trip into the same cells")

cat("connectivity through the raster path\n")

s2 <- adaptels_raster(rr, output = NULL, threshold = 200,
                      queen_topology = TRUE, connectivity = TRUE)
s1 <- adaptels_raster(rr, output = NULL, threshold = 200,
                      queen_topology = TRUE)
ok(s2$n_adaptels >= s1$n_adaptels,
   "connectivity = TRUE never reduces the count at min_size 0")

cat("errors\n")

ok(inherits(try(adaptels_raster(rr, output = NULL, threshold = 60,
                                distance = "nope"), silent = TRUE),
            "try-error"),
   "an unknown distance is rejected through the raster path")
ok(inherits(try(adaptels_raster(mk(nl = 1L), output = NULL, threshold = 60,
                                distance = "cosine"), silent = TRUE),
            "try-error"),
   "an out-of-scale threshold is rejected through the raster path")

cat("\nall raster bridge tests passed\n")
}
