# grow_seeds. Plain-R tests, no testthat dependency.
#
# The point->pixel contract is tested first and on its own (SPEC_grow_seeds
# section 9.9): it is the one place this package and plGeoAdaptels can drift
# silently. The same coordinate table is checked in the Python suite; here the
# expected indices are the Python 0-based (row, col) + 1, because R seeds are
# 1-based (matching sicle()).

library(rgeoadaptels)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

# Same north-up transform as the Python table: EPSG:2180-ish easting/northing,
# 0.25 m pixels, 8 cols x 6 rows. Every coordinate is exactly representable.
X0 <- 500000.0; Y0 <- 400000.0
W <- 0.25; H <- 0.25
COLS <- 8L; ROWS <- 6L

p2p <- rgeoadaptels:::.point_to_pixel

cat("point -> pixel contract (1-based)\n")

# label, x, y, expected_row (1-based), expected_col (1-based)
tbl <- list(
  list("centre of (1,1)",              500000.125, 399999.875, 1L, 1L),
  list("centre of (3,4)",              500000.875, 399999.375, 3L, 4L),
  list("vertical edge -> right pixel", 500001.0,   399999.9,   1L, 5L),
  list("horizontal edge -> pixel below", 500000.1, 399999.5,   3L, 1L),
  list("top-left outer corner -> (1,1)", 500000.0, 400000.0,   1L, 1L),
  list("last valid pixel (6,8) centre", 500001.875, 399998.625, 6L, 8L),
  list("outer east edge -> col ROWS+ (outside)", 500002.0, 399999.9, 1L, 9L),
  list("outer south edge -> row ROWS+ (outside)", 500000.1, 399998.5, 7L, 1L)
)

for (t in tbl) {
  rc <- p2p(t[[2]], t[[3]], X0, Y0, W, H)
  ok(rc[["row"]] == t[[4]] && rc[["col"]] == t[[5]],
     sprintf("%s -> (%d,%d)", t[[1]], t[[4]], t[[5]]))
}

cat("returns integers\n")
rc <- p2p(500000.1, 399999.9, X0, Y0, W, H)
ok(is.integer(rc), "row/col are integer, not double")

cat("edge rule is right-and-below\n")
rc <- p2p(X0 + 3 * W, Y0 - 4.5 * H, X0, Y0, W, H)
ok(rc[["col"]] == 4L, "east edge belongs to the right pixel (col 4, 1-based)")
rc <- p2p(X0 + 1.5 * W, Y0 - 3 * H, X0, Y0, W, H)
ok(rc[["row"]] == 4L, "south edge belongs to the pixel below (row 4, 1-based)")

cat("brute-force floor over every pixel centre\n")
allok <- TRUE
for (r in seq_len(ROWS)) for (cc in seq_len(COLS)) {
  x <- X0 + (cc - 1 + 0.5) * W
  y <- Y0 - (r - 1 + 0.5) * H
  rc <- p2p(x, y, X0, Y0, W, H)
  if (rc[["row"]] != r || rc[["col"]] != cc) allok <- FALSE
}
ok(allok, "every pixel centre round-trips to its own (row, col)")

# ── growth (SPEC 9.1-9.5 + the section 3 features), mirroring the Python ──
# suite. Seeds here are 1-based; the Python (r, c) become (r+1, c+1).

two_blocks <- function(rows = 12L, cols = 20L, mid = 10L, lo = 10, hi = 200) {
  d <- array(0, dim = c(1L, rows, cols))
  d[1, , seq_len(mid)] <- lo
  d[1, , (mid + 1L):cols] <- hi
  d
}

cat("two uniform blocks, one seed each\n")
d <- two_blocks()
o <- grow_seeds(d, rbind(c(4, 1), c(4, 20)), quiet = TRUE)
ok(identical(sort(unique(as.vector(o$labels))), c(-1L, 0L, 1L)) == FALSE &&
   identical(sort(unique(as.vector(o$labels))), c(0L, 1L)),
   "two segments, nothing unassigned")
ok(all(o$labels[, 1:10] == 0L) && all(o$labels[, 11:20] == 1L),
   "boundary exactly on the block edge")

cat("max_cost confines a single seed to its block\n")
o <- grow_seeds(two_blocks(), rbind(c(4, 1)), max_cost = 50, quiet = TRUE)
ok(all(o$labels[, 1:10] == 0L) && all(o$labels[, 11:20] == -1L),
   "its block only; the far block stays -1")

cat("compactness splits a uniform image at the midpoint\n")
u <- array(42, dim = c(1L, 15L, 20L))
o <- grow_seeds(u, rbind(c(8, 5), c(8, 16)), compactness = 1, quiet = TRUE)
ok(all(o$labels[, 1:10] == 0L) && all(o$labels[, 11:20] == 1L),
   "boundary at the geometric midpoint, deterministically")

cat("label i is seeds[i, ] under shuffle\n")
set.seed(4)
d3 <- array(runif(3 * 30 * 30, 0, 255), dim = c(3, 30, 30))
base_seeds <- rbind(c(6, 6), c(6, 26), c(26, 6), c(26, 26), c(16, 16))
allok <- TRUE
for (rep in 1:4) {
  s <- base_seeds[sample(nrow(base_seeds)), ]
  o <- grow_seeds(d3, s, quiet = TRUE)
  for (i in seq_len(nrow(s)))
    if (o$labels[s[i, 1], s[i, 2]] != (i - 1L)) allok <- FALSE
}
ok(allok, "each seed pixel carries its own 0-based index, any order")

cat("defaults leave nothing unassigned (kernel partition)\n")
o <- grow_seeds(d3, base_seeds, quiet = TRUE)
ok(sum(o$labels == -1L) == 0L, "no cap, no radius -> full partition")

cat("band_weights can switch a band off\n")
db <- array(0, dim = c(2L, 12L, 20L))
db[1, , ] <- 100
db[2, , 1:10] <- 0; db[2, , 11:20] <- 100
on  <- grow_seeds(db, rbind(c(7, 1)), max_cost = 50, quiet = TRUE)
off <- grow_seeds(db, rbind(c(7, 1)), max_cost = 50, band_weights = c(1, 0),
                  quiet = TRUE)
ok(all(on$labels[, 11:20] == -1L), "weighted: the step blocks the far half")
ok(all(off$labels[, 11:20] == 0L), "band zeroed: the far half grows in")

cat("seed_window median rescues an outlier click\n")
uo <- array(50, dim = c(1L, 15L, 15L))
uo[1, 8, 8] <- 200
raw <- grow_seeds(uo, rbind(c(8, 8)), max_cost = 50, seed_window = 1, quiet = TRUE)
med <- grow_seeds(uo, rbind(c(8, 8)), max_cost = 50, seed_window = 3, quiet = TRUE)
ok(sum(raw$labels >= 0L) == 1L, "raw pixel anchors on 200; nothing grows")
ok(sum(med$labels >= 0L) == 15L * 15L, "3x3 median = 50; the region grows")

cat("seed_window must be odd\n")
ok(inherits(try(grow_seeds(array(10, c(1L, 8L, 8L)), rbind(c(4, 4)),
                           seed_window = 2, quiet = TRUE), silent = TRUE),
            "try-error"), "even seed_window is refused")

cat("max_radius caps the reach\n")
ur <- array(42, dim = c(1L, 21L, 21L))
o <- grow_seeds(ur, rbind(c(11, 11)), max_radius = 5, quiet = TRUE)
ok(o$labels[11, 14] == 0L, "3 px away is within the radius")
ok(o$labels[11, 19] == -1L, "8 px away is beyond it")
far_ok <- TRUE
for (r in seq_len(21)) for (cc in seq_len(21)) {
  dist <- sqrt((r - 11)^2 + (cc - 11)^2)
  if (dist > 5 && o$labels[r, cc] != -1L) far_ok <- FALSE
}
ok(far_ok, "everything past the radius is unassigned")

cat("fill_holes closes an interior pocket\n")
uh <- array(50, dim = c(1L, 11L, 11L)); uh[1, 4, 4] <- 200   # interior highlight
open_ <- grow_seeds(uh, rbind(c(6, 6)), max_cost = 30, quiet = TRUE)
ok(open_$labels[4, 4] == -1L, "highlight over the cap -> a hole")
filled <- grow_seeds(uh, rbind(c(6, 6)), max_cost = 30, fill_holes = TRUE,
                     quiet = TRUE)
ok(filled$labels[4, 4] == 0L, "fill_holes closes the enclosed pocket")
ok(sum(filled$labels != open_$labels) == 1L, "only the pocket changed")

cat("fill_holes leaves nodata and edges alone\n")
un <- array(50, dim = c(1L, 11L, 11L)); un[1, 4, 5] <- 200
mk <- matrix(0L, 11L, 11L); mk[4, 4] <- 1L                   # nodata beside it
fn <- grow_seeds(un, rbind(c(6, 6)), mask = mk, max_cost = 30,
                 fill_holes = TRUE, quiet = TRUE)
ok(fn$labels[4, 5] == -1L, "a pocket touching nodata is not filled")
ok(fn$labels[4, 4] == -1L, "nodata itself stays unassigned")

# ── file wrapper (SPEC 4.2, 5, 9.6): point layer -> pixels -> growth on a ──
# terra raster, with a real vector round-trip and a CRS reprojection. Skipped
# where terra is absent, like the rest of the raster path.
if (requireNamespace("terra", quietly = TRUE)) {
  cat("grow_seeds_raster: point layer on a terra raster\n")
  set.seed(1)
  nr <- 40L; nc <- 50L; rs <- 0.25
  r <- terra::rast(nrows = nr, ncols = nc, xmin = 500000,
                   xmax = 500000 + nc * rs, ymin = 400000 - nr * rs,
                   ymax = 400000, crs = "EPSG:2180", nlyrs = 3)
  terra::values(r) <- matrix(sample(0:255, nr * nc * 3, TRUE), ncol = 3)
  x0 <- terra::xmin(r); y0 <- terra::ymax(r)
  w <- terra::xres(r); h <- terra::yres(r)
  px <- rbind(c(10, 10), c(10, 40), c(30, 10), c(30, 40), c(20, 25))
  xy <- t(apply(px, 1, function(p) c(x0 + (p[2] - 0.5) * w,
                                     y0 - (p[1] - 0.5) * h)))
  v <- terra::vect(xy, type = "points", crs = "EPSG:2180")

  seg <- grow_seeds_raster(r, v, quiet = TRUE)
  ok(all(vapply(seq_len(nrow(px)), function(i)
    seg$labels[px[i, 1], px[i, 2]] == (i - 1L), logical(1))),
    "each point lands on its own pixel (label i-1)")

  v84 <- terra::project(v, "EPSG:4326")
  seg2 <- grow_seeds_raster(r, v84, quiet = TRUE)
  ok(all(vapply(seq_len(nrow(px)), function(i)
    seg2$labels[px[i, 1], px[i, 2]] == (i - 1L), logical(1))),
    "points in EPSG:4326 reproject to the same pixels")

  td <- tempfile(); dir.create(td)
  ot <- file.path(td, "labels.tif"); op <- file.path(td, "crowns.gpkg")
  grow_seeds_raster(r, v, output = ot, polygons = op, max_cost = 40,
                    quiet = TRUE)
  ok(file.exists(ot) && file.exists(op), "label raster and polygons written")
  p <- terra::vect(op)
  ids <- terra::values(p)[[1]]
  ok(nrow(p) >= 1L && all(ids >= 0 & ids < nrow(px)),
     "polygons carry valid segment ids")

  cat("matrix of map coordinates works too\n")
  seg3 <- grow_seeds_raster(r, xy, quiet = TRUE)
  ok(identical(seg3$labels, seg$labels),
     "a plain (x, y) matrix matches the SpatVector")
} else {
  cat("terra not installed -- skipping grow_seeds_raster tests\n")
}

cat("ALL GROW TESTS PASSED\n")
