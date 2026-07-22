# SICLE. Plain-R tests, no testthat dependency.
#
# These do not compare against plGeoAdaptels -- that needs Python and lives
# in tools/. They pin what this package promises on its own.

library(rgeoadaptels)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

scene <- function(bands = 3L, rows = 45L, cols = 58L, seed = 2L) {
  set.seed(seed)
  array(as.numeric(sample(0:255, bands * rows * cols, TRUE)),
        dim = c(bands, rows, cols))
}

d <- scene()
ROWS <- 45L; COLS <- 58L

seeds_of <- function(n, seed = 3L) {
  set.seed(seed)
  i <- sample.int(ROWS * COLS, n)
  cbind(((i - 1L) %/% COLS) + 1L, ((i - 1L) %% COLS) + 1L)
}

# ── structure ────────────────────────────────────────────────────────
cat("structure\n")

set.seed(1)
o <- sicle(d, n_segments = 40, quiet = TRUE)
ok(identical(dim(o$labels), c(ROWS, COLS)), "labels keep the raster's shape")
ok(is.integer(o$labels), "labels are integer")
ok(o$n_superpixels == 40L, "n_segments is met exactly")
ok(identical(sort(unique(as.vector(o$labels))), 0:39),
   "ids are consecutive from 0")
ok(all(o$labels >= 0), "every pixel is assigned")
ok(nrow(o$labels) != ncol(o$labels), "the test scene is not square")

for (n in c(5L, 40L, 150L)) {
  set.seed(1)
  ok(sicle(d, n_segments = n, quiet = TRUE)$n_superpixels == n,
     paste("n_segments =", n, "is met exactly"))
}

set.seed(7)
o2 <- sicle(matrix(as.numeric(sample(0:255, 30 * 37, TRUE)), 30, 37),
            n_segments = 12, quiet = TRUE)
ok(identical(dim(o2$labels), c(30L, 37L)), "a 2-D matrix is accepted")

# ── the iteration count does less than it looks ──────────────────────
cat("iterations\n")

s <- seeds_of(700)
a2 <- sicle(d, seeds = s, n_segments = 60, n_iterations = 2, quiet = TRUE)
a3 <- sicle(d, seeds = s, n_segments = 60, n_iterations = 3, quiet = TRUE)
a5 <- sicle(d, seeds = s, n_segments = 60, n_iterations = 5, quiet = TRUE)

# The paper's curve is M(i) = max(N0^(1 - i/(Omega-1)), Nf), whose exponent
# never sees n_segments. With N0 = 700 and Nf = 60 that makes Omega = 3
# perform exactly one removal step, the same as Omega = 2. Documented here
# because it looks like a bug and is not.
ok(identical(a2$labels, a3$labels),
   "n_iterations = 3 is bit-identical to 2, as the paper's curve implies")
ok(!identical(a2$labels, a5$labels), "n_iterations = 5 does differ")
ok(a2$n_superpixels == 60L && a5$n_superpixels == 60L,
   "the target is met whatever the iteration count")

# ── seeds ────────────────────────────────────────────────────────────
cat("seeds\n")

ok(sicle(d, seeds = s, n_segments = 60, quiet = TRUE)$n_superpixels == 60L,
   "explicit seeds are accepted")

b1 <- sicle(d, seeds = s, n_segments = 60, quiet = TRUE)
b2 <- sicle(d, seeds = s, n_segments = 60, quiet = TRUE)
ok(identical(b1$labels, b2$labels), "explicit seeds make the run reproducible")

# n_oversampling = 500 on purpose. The scene holds 2610 pixels, so the
# default of 3000 makes every pixel a seed: the "sample" is the whole raster
# and only its order changes, which the result does not depend on. A test of
# set.seed() written on the default would pass for the wrong reason.
set.seed(11); c1 <- sicle(d, n_segments = 40, n_oversampling = 500, quiet = TRUE)
set.seed(11); c2 <- sicle(d, n_segments = 40, n_oversampling = 500, quiet = TRUE)
ok(identical(c1$labels, c2$labels), "set.seed() makes sampling reproducible")

set.seed(12); c3 <- sicle(d, n_segments = 40, n_oversampling = 500, quiet = TRUE)
ok(!identical(c1$labels, c3$labels), "a different seed gives a different result")

set.seed(11); e1 <- sicle(d, n_segments = 40, quiet = TRUE)
set.seed(12); e2 <- sicle(d, n_segments = 40, quiet = TRUE)
ok(identical(e1$labels, e2$labels),
   "the seed is inert once N0 covers every valid pixel")

err <- function(expr, what) {
  ok(inherits(try(expr, silent = TRUE), "try-error"), what)
}
err(sicle(d, seeds = cbind(1, 1, 1), n_segments = 2, quiet = TRUE),
    "a three-column seed matrix is rejected")
err(sicle(d, seeds = rbind(c(1, 1), c(999, 1)), n_segments = 2, quiet = TRUE),
    "a seed outside the raster is rejected")
err(sicle(d, seeds = rbind(c(5, 5), c(5, 5), c(9, 9)), n_segments = 2,
          quiet = TRUE),
    "duplicate seeds are rejected")
err(sicle(d, seeds = s[1:10, ], n_segments = 60, quiet = TRUE),
    "fewer seeds than segments is rejected")

m <- matrix(0L, ROWS, COLS); m[1:5, ] <- 1L
err(sicle(d, mask = m, seeds = rbind(c(1, 1), c(2, 2), c(30, 30)),
          n_segments = 2, quiet = TRUE),
    "seeds on nodata are rejected")

# ── nodata ───────────────────────────────────────────────────────────
cat("nodata\n")

set.seed(5)
om <- sicle(d, mask = m, n_segments = 30, quiet = TRUE)
ok(all(om$labels[m == 1] == -1L), "masked pixels come back as -1")
ok(all(om$labels[m == 0] >= 0), "unmasked pixels are all assigned")
ok(identical(sort(unique(as.vector(om$labels[m == 0]))), 0:29),
   "ids stay consecutive when a mask is present")

dna <- d; dna[1, 6, 6] <- NA
set.seed(5)
ok(sicle(dna, n_segments = 20, quiet = TRUE)$labels[6, 6] == -1L,
   "NA in a band is treated as nodata when no mask is given")

# ── saliency ─────────────────────────────────────────────────────────
cat("saliency\n")

sal <- matrix(0, ROWS, COLS); sal[10:25, 15:35] <- 1
with_sal <- sicle(d, seeds = s, n_segments = 60, saliency = sal, quiet = TRUE)
without <- sicle(d, seeds = s, n_segments = 60, quiet = TRUE)
ok(with_sal$n_superpixels == 60L, "saliency does not change the target")
ok(!identical(with_sal$labels, without$labels),
   "saliency actually steers seed removal")

sal_na <- sal; sal_na[3, 3] <- NA
err(sicle(d, seeds = s, n_segments = 60, saliency = sal_na, quiet = TRUE),
    "NA saliency over a valid pixel is rejected")

# NA under nodata is fine: it never reaches the relevance.
sal_ok <- sal; sal_ok[1:5, ] <- NA
s_valid <- s[m[s] == 0L | TRUE, , drop = FALSE]
s_valid <- s_valid[m[cbind(s_valid[, 1], s_valid[, 2])] == 0L, , drop = FALSE]
ok(sicle(d, mask = m, seeds = s_valid, n_segments = 30, saliency = sal_ok,
         quiet = TRUE)$n_superpixels == 30L,
   "NA saliency under nodata is accepted")

err(sicle(d, seeds = s, n_segments = 60, saliency = matrix(0, 3, 3),
          quiet = TRUE),
    "a saliency map of the wrong size is rejected")

# ── connectivity of the result ───────────────────────────────────────
cat("connectivity\n")

# SICLE grows over an 8-adjacency, so every label must be one 8-connected
# region. Deliberately not enforce_connectivity(), which tests 4 and would
# report a defect that is not there -- the same trap documented for the
# Python package.
lab <- b1$labels
eight_connected <- function(lab, id) {
  m <- lab == id
  idx <- which(m, arr.ind = TRUE)
  seen <- matrix(FALSE, nrow(lab), ncol(lab))
  stack <- idx[1, , drop = FALSE]
  seen[stack[1, 1], stack[1, 2]] <- TRUE
  n <- 1L
  while (nrow(stack)) {
    p <- stack[1, ]; stack <- stack[-1, , drop = FALSE]
    for (dr in -1:1) for (dc in -1:1) {
      r <- p[1] + dr; c <- p[2] + dc
      if (r < 1 || r > nrow(lab) || c < 1 || c > ncol(lab)) next
      if (!m[r, c] || seen[r, c]) next
      seen[r, c] <- TRUE; n <- n + 1L
      stack <- rbind(stack, c(r, c))
    }
  }
  n == sum(m)
}
split <- Filter(function(id) !eight_connected(lab, id),
                sort(unique(as.vector(lab))))
ok(length(split) == 0L, "every superpixel is a single 8-connected region")

# ── validation ───────────────────────────────────────────────────────
cat("validation\n")

err(sicle("not a raster", n_segments = 10, quiet = TRUE),
    "non-numeric data is rejected")
err(sicle(array(0, dim = c(2, 2, 2, 2)), n_segments = 2, quiet = TRUE),
    "a 4-D array is rejected")
err(sicle(d, n_segments = 0, quiet = TRUE), "n_segments = 0 is rejected")
err(sicle(d, mask = matrix(0, 3, 3), n_segments = 10, quiet = TRUE),
    "a mask of the wrong size is rejected")
err(sicle(d, n_segments = ROWS * COLS + 1L, quiet = TRUE),
    "asking for more segments than pixels is rejected")

got_warning <- FALSE
withCallingHandlers(
  invisible(sicle(d, n_segments = 30, n_oversampling = 5, quiet = TRUE)),
  warning = function(w) {
    got_warning <<- grepl("n_oversampling", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
ok(got_warning, "n_oversampling below n_segments warns rather than correcting quietly")

cat("\nall SICLE tests passed\n")
