# grow_seeds -- seeded spectral region growing ("inverse OBIA").
#
# Grow each operator-placed point into the region that looks like the pixel it
# sits on, leaving everything unseeded unassigned (-1). The growth is the SICLE
# IFT core, unchanged: C_ift_fmax is called once, with all seeds, and every
# option is prepared into its inputs or read out of its outputs. The kernel
# stays frozen so this twin and the Python one share one verified algorithm and
# differ only in this wrapper -- see docs/SPEC_grow_seeds.md section 3.
#
# Every step here mirrors plgeoadaptels/grow.py operation for operation and in
# the same order (weights, then the seed-window median on the weighted buffer,
# then the coordinate bands), because the two must agree bit for bit -- that is
# what tools/cross_validate_against_plgeoadaptels.R checks.
#
# Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.


# Map a map coordinate to a 1-based (row, col) pixel index.
#
# Contract (SPEC_grow_seeds section 5.1): pixel (r, c) covers the half-open
# extent [x0 + c*w, x0 + (c+1)*w) x (y0 - r*h, y0 - (r+1)*h), so a point on a
# shared grid line belongs to the pixel to the right / below. floor gives
# exactly that. Written out rather than delegated to terra::cellFromXY, so the
# specification is the contract and this matches Python's _point_to_pixel bit
# for bit (SPEC 5.2). Python is 0-based; this returns the same (r, c) + 1 to
# match sicle()'s (row, col) seed convention. x0/y0 are the top-left corner;
# w, h are positive pixel size. May fall outside the raster; bounds are checked
# separately (SPEC 5.4).
.point_to_pixel <- function(x, y, x0, y0, w, h) {
  col0 <- floor((x - x0) / w)
  row0 <- floor((y0 - y) / h)
  c(row = as.integer(row0) + 1L, col = as.integer(col0) + 1L)
}


# 4-connected dilation of a logical matrix, one step, confined to `within`.
# Vectorised matrix shifts -- no per-pixel loop -- so it is fast enough to
# iterate to convergence on a full tile.
.dilate4 <- function(m, within, rows, cols) {
  d <- m
  if (rows > 1L) {
    d[2:rows, ] <- d[2:rows, ] | m[1:(rows - 1L), ]
    d[1:(rows - 1L), ] <- d[1:(rows - 1L), ] | m[2:rows, ]
  }
  if (cols > 1L) {
    d[, 2:cols] <- d[, 2:cols] | m[, 1:(cols - 1L)]
    d[, 1:(cols - 1L)] <- d[, 1:(cols - 1L)] | m[, 2:cols]
  }
  d & within
}

.flood <- function(seed_mask, within, rows, cols) {
  keep <- seed_mask & within
  repeat {
    nk <- .dilate4(keep, within, rows, cols)
    if (all(nk == keep)) break
    keep <- nk
  }
  keep
}

# Fill unassigned pockets that sit fully inside one crown: a 4-connected group
# of valid -1 pixels that does not reach the border, does not touch nodata,
# and is bounded by exactly one label. Conservative on purpose -- those
# conditions keep it from swallowing the nodata region or bridging crowns --
# and topological, so the Python twin reaches the same result.
.fill_holes <- function(labels, mask, rows, cols) {
  out <- labels
  fillable <- (labels < 0L) & (mask == 0L)
  if (!any(fillable)) return(out)
  border <- matrix(FALSE, rows, cols)
  border[1, ] <- TRUE; border[rows, ] <- TRUE
  border[, 1] <- TRUE; border[, cols] <- TRUE
  reach <- .flood(border & fillable, fillable, rows, cols)
  todo <- fillable & !reach
  while (any(todo)) {
    idx <- which(todo)[1]
    r0 <- ((idx - 1L) %% rows) + 1L; c0 <- ((idx - 1L) %/% rows) + 1L
    seedm <- matrix(FALSE, rows, cols); seedm[r0, c0] <- TRUE
    comp <- .flood(seedm, todo, rows, cols)
    nb <- .dilate4(comp, matrix(TRUE, rows, cols), rows, cols) & !comp
    if (!any(mask[nb] != 0L)) {
      nl <- unique(labels[nb & (labels >= 0L)])
      if (length(nl) == 1L) out[comp] <- nl
    }
    todo <- todo & !comp
  }
  out
}


#' Seeded spectral region growing ("inverse OBIA")
#'
#' Grow each operator-placed point into the region that looks like the pixel it
#' sits on, leaving everything unseeded unassigned (`-1`). This is the inverse
#' of [adaptels()] and [sicle()], which partition the whole image: here the
#' operator supplies the objects and the algorithm supplies their boundaries.
#'
#' The growth is the SICLE IFT core called once, with every seed, and no seed
#' is ever removed, so `labels == i` is exactly the region grown from the
#' `i`-th seed and `-1` is unassigned. This is the R twin of plGeoAdaptels'
#' `grow_seeds`; the two produce identical labels on identical input.
#'
#' @param data Numeric matrix `(rows, cols)` for one band, or array
#'   `(bands, rows, cols)` for several. Bands are treated as a feature vector;
#'   what they are (RGB, CIELAB, CIR, an index, ...) is the operator's choice,
#'   made upstream.
#' @param seeds `(n, 2)` matrix of `(row, col)` starting points, **1-based**,
#'   matching [sicle()]. One row per object of interest.
#' @param mask Optional matrix `(rows, cols)`. Non-zero marks nodata; `NA` in
#'   any band is treated as nodata when `mask` is `NULL`.
#' @param max_cost Cost cap in band units -- a tolerance on the minimax
#'   spectral deviation from the seed. `NULL` keeps every reachable pixel. If
#'   the input was CIELAB this is a Delta-E tolerance.
#' @param band_weights Optional length-`bands` multipliers applied before the
#'   distance.
#' @param compactness SLIC-style spatial term in feature-units per pixel; `0`
#'   reproduces pure spectral growth exactly.
#' @param seed_window Odd side `k` of a `k*k` median used as the seed
#'   signature; `1` is the raw pixel.
#' @param max_radius Hard limit in pixels; pixels further than this from their
#'   seed go back to `-1`.
#' @param fill_holes Fill unassigned pockets that sit fully inside one crown --
#'   the interior pixels a `max_cost` cut left as `-1`. A pocket is filled only
#'   when it does not reach the border, does not touch nodata, and is bounded
#'   by a single label, so it cannot swallow the nodata region or bridge two
#'   crowns.
#' @param return_cost Return the per-pixel path cost as well.
#' @param quiet Suppress the progress message.
#'
#' @return A list with `labels`, an integer matrix of 0-based labels (`-1` for
#'   unassigned) where label `i` is the region grown from `seeds[i, ]`;
#'   `n_segments`, the number of seeds; and `cost`, a numeric matrix when
#'   `return_cost` is `TRUE`, else `NULL`.
#'
#' @seealso [sicle()] and [adaptels()] for unsupervised segmentation.
#' @examples
#' d <- array(runif(3 * 40 * 40, 0, 255), dim = c(3, 40, 40))
#' seeds <- rbind(c(10, 10), c(30, 30))
#' out <- grow_seeds(d, seeds, max_cost = 60, quiet = TRUE)
#' out$n_segments
#' @export
grow_seeds <- function(data, seeds, mask = NULL, max_cost = NULL,
                       band_weights = NULL, compactness = 0,
                       seed_window = 1L, max_radius = NULL,
                       fill_holes = FALSE,
                       return_cost = FALSE, quiet = FALSE) {
  if (!is.numeric(data)) stop("data must be numeric", call. = FALSE)
  d <- dim(data)
  if (length(d) == 2L) {
    n_layers <- 1L; rows <- d[1]; cols <- d[2]
    arr <- array(data, dim = c(1L, rows, cols))
  } else if (length(d) == 3L) {
    n_layers <- d[1]; rows <- d[2]; cols <- d[3]
    arr <- data
  } else {
    stop("data must be a 2-D matrix or a 3-D (bands, rows, cols) array",
         call. = FALSE)
  }
  cols <- as.integer(cols); rows <- as.integer(rows)
  size <- rows * cols
  # Row-major, band-sequential, matching the NumPy (n_layers, size) C-order
  # buffer -- flat[l*size + r*cols + c]. Same flatten sicle() uses.
  flat <- as.double(as.vector(aperm(arr, c(3, 2, 1))))

  if (is.null(mask)) {
    m <- rep(0L, size)
    for (l in seq_len(n_layers))
      m <- m | as.integer(is.na(flat[(l - 1) * size + seq_len(size)]))
    m <- as.integer(m)
  } else {
    if (!identical(dim(as.matrix(mask)), c(rows, cols)))
      stop("mask must be ", rows, "x", cols, call. = FALSE)
    m <- as.integer(as.vector(t(as.matrix(mask))) != 0)
  }
  flat[is.na(flat)] <- 0    # masked values are never read; zero for parity

  # Seeds: 1-based (row, col) -> 0-based flat indices, validated. Raise rather
  # than drop -- dropping would renumber labels and break the point-order
  # contract (SPEC 4.3), exactly as the Python twin does.
  seeds <- as.matrix(seeds)
  if (ncol(seeds) != 2L)
    stop("seeds must be an (n, 2) matrix of (row, col) pairs, got ",
         ncol(seeds), " column(s)", call. = FALSE)
  if (nrow(seeds) == 0L)
    stop("seeds is empty; grow_seeds needs at least one point", call. = FALSE)
  sr <- as.integer(seeds[, 1]); sc <- as.integer(seeds[, 2])
  outside <- sr < 1L | sr > rows | sc < 1L | sc > cols
  if (any(outside))
    stop("seed(s) at input index ", paste(which(outside), collapse = ", "),
         " lie outside the raster (", rows, "x", cols, ")", call. = FALSE)
  seed_flat <- (sr - 1L) * cols + (sc - 1L)          # 0-based
  on_nodata <- m[seed_flat + 1L] != 0L
  if (any(on_nodata))
    stop("seed(s) at input index ", paste(which(on_nodata), collapse = ", "),
         " fall on nodata pixels", call. = FALSE)
  if (anyDuplicated(seed_flat)) {
    dup <- which(duplicated(seed_flat) | duplicated(seed_flat, fromLast = TRUE))
    stop("seed(s) at input index ", paste(dup, collapse = ", "),
         " land on the same pixel(s); two points in one pixel would merge two ",
         "segments and break the label contract", call. = FALSE)
  }

  if (!is.null(band_weights)) {
    bw <- as.double(band_weights)
    if (length(bw) != n_layers)
      stop("band_weights must have one entry per band (", n_layers, "), got ",
           length(bw), call. = FALSE)
    for (l in seq_len(n_layers)) {
      sl <- (l - 1) * size + seq_len(size)
      flat[sl] <- flat[sl] * bw[l]
    }
  }

  seed_window <- as.integer(seed_window)
  if (seed_window > 1L) {
    if (seed_window %% 2L == 0L)
      stop("seed_window must be odd, got ", seed_window,
           ": an even window has no centre pixel to anchor on", call. = FALSE)
    half <- seed_window %/% 2L
    # All medians from the pre-write buffer, then written, so one seed's window
    # cannot be perturbed by another seed's overwrite. Median on the weighted
    # buffer, matching the units the kernel compares against.
    meds <- vector("list", length(seed_flat))
    for (si in seq_along(seed_flat)) {
      pr <- sr[si] - 1L; pc <- sc[si] - 1L            # 0-based
      r0 <- max(0L, pr - half); r1 <- min(rows - 1L, pr + half)
      c0 <- max(0L, pc - half); c1 <- min(cols - 1L, pc + half)
      wr <- rep(r0:r1, each = (c1 - c0 + 1L))
      wc <- rep(c0:c1, times = (r1 - r0 + 1L))
      widx <- wr * cols + wc                          # 0-based flat
      widx <- widx[m[widx + 1L] == 0L]                # unmasked only
      med <- numeric(n_layers)
      for (l in seq_len(n_layers))
        med[l] <- stats::median(flat[(l - 1) * size + widx + 1L])
      meds[[si]] <- med
    }
    for (si in seq_along(seed_flat))
      for (l in seq_len(n_layers))
        flat[(l - 1) * size + seed_flat[si] + 1L] <- meds[[si]][l]
  }

  n_used <- n_layers
  if (compactness > 0) {
    # Two coordinate bands scaled by lambda; because wroot is seed-anchored the
    # appended terms evaluate to lambda^2 * d_euclid(seed, pixel)^2 inside the
    # norm -- SLIC-style compactness with no kernel change (SPEC 3.3). Same
    # row-major order and values as the Python np.mgrid ravel.
    idx0 <- seq_len(size) - 1L
    rr <- idx0 %/% cols
    ccx <- idx0 %% cols
    flat <- c(flat, compactness * rr, compactness * ccx)
    n_used <- n_layers + 2L
  }

  res <- .Call(C_ift_fmax, flat, m, cols, rows, seed_flat)
  labels <- res[[1]]; cost <- res[[2]]

  if (!is.null(max_cost)) labels[cost > max_cost] <- -1L

  if (!is.null(max_radius)) {
    assigned <- which(labels >= 0L)
    ai <- assigned - 1L                               # 0-based flat
    pr <- ai %/% cols; pc <- ai %% cols
    lab <- labels[assigned]                           # 0-based seed index
    dr <- pr - (sr[lab + 1L] - 1L)
    dc <- pc - (sc[lab + 1L] - 1L)
    dist <- sqrt(dr * dr + dc * dc)
    labels[assigned[dist > max_radius]] <- -1L
  }

  lab_mat <- matrix(labels, nrow = rows, ncol = cols, byrow = TRUE)

  # Fill unassigned pockets left inside a crown by the max_cost cut, after the
  # cut that creates them.
  if (isTRUE(fill_holes)) {
    mask_mat <- matrix(m, nrow = rows, ncol = cols, byrow = TRUE)
    lab_mat <- .fill_holes(lab_mat, mask_mat, rows, cols)
  }

  if (!quiet) {
    n_assigned <- sum(lab_mat >= 0L)
    message(sprintf("grow_seeds: %d seeds, %d px assigned, %d px unassigned",
                    length(seed_flat), n_assigned, size - n_assigned))
  }

  list(labels = lab_mat,
       n_segments = length(seed_flat),
       cost = if (return_cost)
         matrix(cost, nrow = rows, ncol = cols, byrow = TRUE) else NULL)
}


#' Grow a point layer into regions on a raster from disk
#'
#' Reads the raster and the point layer, converts each point to the pixel it
#' falls on (SPEC_grow_seeds section 5), grows every point in one call, and
#' writes a label raster and/or crown polygons. `...` is passed to
#' [grow_seeds()] (`max_cost`, `band_weights`, `compactness`, `seed_window`,
#' `max_radius`).
#'
#' @param input A `SpatRaster`, a path, or several paths taken as bands. Feed
#'   the bands the growth should see (RGB, CIELAB, an index), prepared upstream.
#' @param points A `SpatVector` of points, a path to any OGR-readable point
#'   layer, or a two-column matrix of `(x, y)` map coordinates already in the
#'   raster CRS. A layer whose CRS differs from the raster's is reprojected,
#'   and that is reported unless `quiet` -- never assumed.
#' @param output Path for the `int32` label raster (`-1` = unassigned, set as
#'   its NA value). `NULL` writes none.
#' @param polygons Path for crown polygons; the format follows the extension
#'   (`.gpkg` recommended -- one file, no field-name limits, so the join back
#'   to point attributes is clean). `NULL` writes none.
#' @param quiet Suppress progress messages.
#' @param ... Passed to [grow_seeds()].
#'
#' @return Invisibly, the list [grow_seeds()] returns. Label `i` is the region
#'   grown from the i-th point, in feature order, so it joins to that point.
#'
#' @details Requires \pkg{terra}, which reads the raster and the vector and
#'   does the reprojection. The growth itself has no spatial dependency.
#'
#' @seealso [grow_seeds()] for the array interface, [sicle()] and [adaptels()]
#'   for unsupervised segmentation.
#' @export
grow_seeds_raster <- function(input, points, output = NULL, polygons = NULL,
                              quiet = FALSE, ...) {
  .need_terra()
  b <- read_bands(input)
  tmpl <- b$template
  x0 <- terra::xmin(tmpl); y0 <- terra::ymax(tmpl)
  w <- terra::xres(tmpl);  h <- terra::yres(tmpl)

  if (is.matrix(points) || is.data.frame(points)) {
    xy <- as.matrix(points)
    if (ncol(xy) != 2L)
      stop("points matrix must be (n, 2) of (x, y) map coordinates in the ",
           "raster CRS", call. = FALSE)
  } else {
    v <- if (inherits(points, "SpatVector")) points else terra::vect(points)
    if (terra::geomtype(v) != "points")
      stop("grow_seeds needs a point layer; got ", terra::geomtype(v),
           ". Digitise the objects as points.", call. = FALSE)
    vc <- terra::crs(v); tc <- terra::crs(tmpl)
    if (nzchar(vc) && nzchar(tc) && !terra::same.crs(v, tmpl)) {
      v <- terra::project(v, tmpl)
      if (!quiet)
        message(sprintf("  Reprojected %d points to the raster CRS",
                        nrow(terra::crds(v))))
    }
    xy <- terra::crds(v)
  }

  # Each row -> 1-based (row, col). apply returns 2 x n; transpose to n x 2.
  seeds <- t(apply(xy, 1L, function(p)
    .point_to_pixel(p[1], p[2], x0, y0, w, h)))
  seeds <- matrix(as.integer(seeds), ncol = 2L,
                  dimnames = list(NULL, c("row", "col")))

  seg <- grow_seeds(b$data, seeds, mask = b$mask, quiet = quiet, ...)

  if (!is.null(output)) {
    r <- terra::rast(tmpl[[1]])
    # t(): terra fills a layer row-major and labels is already row-major.
    terra::values(r) <- as.vector(t(seg$labels))
    names(r) <- "segment"
    terra::NAflag(r) <- -1
    terra::writeRaster(r, output, overwrite = TRUE, datatype = "INT4S")
    if (!quiet) message(sprintf("  Wrote labels: %s", output))
  }

  if (!is.null(polygons)) {
    r <- terra::rast(tmpl[[1]])
    lab <- seg$labels
    lab[lab < 0L] <- NA                       # only grown regions become polys
    terra::values(r) <- as.vector(t(lab))
    names(r) <- "segment"
    p <- terra::as.polygons(r, dissolve = TRUE, na.rm = TRUE)
    terra::writeVector(p, polygons, overwrite = TRUE)
    if (!quiet)
      message(sprintf("  Wrote %d polygons: %s", nrow(p), polygons))
  }

  invisible(seg)
}
