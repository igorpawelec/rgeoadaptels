#' SICLE superpixels
#'
#' Superpixels through Iterative CLEarcutting: start from far more seeds
#' than you want, grow an optimum-path forest, score every seed, discard the
#' least relevant, and repeat. Unlike [adaptels()] the count is a target you
#' set, not an outcome of the scene.
#'
#' @param data Numeric matrix `(rows, cols)` for one band, or array
#'   `(bands, rows, cols)` for several.
#' @param mask Optional matrix `(rows, cols)`. Non-zero marks nodata. When
#'   `NULL`, `NA` in any band is treated as nodata.
#' @param n_segments Superpixels wanted.
#' @param n_oversampling Seeds to start from. Ignored when `seeds` is given.
#' @param n_iterations Maximum iterations, the paper's `Omega`. See the note
#'   below before trusting the default of 2.
#' @param saliency Optional matrix `(rows, cols)` in [0, 1]. Seeds near
#'   saliency borders survive removal. Must not be `NA` where the mask calls
#'   the pixel valid.
#' @param seeds Optional `(n, 2)` matrix of starting seeds as (row, col),
#'   **1-based**. When given, `n_oversampling` is unused.
#' @param quiet Suppress progress messages.
#'
#' @return A list with `labels`, an integer matrix of 0-based superpixel ids
#'   with `-1` for nodata, and `n_superpixels`.
#'
#' @details
#' `n_iterations` does less than it looks. The paper's preservation curve is
#' `M(i) = max(N0^(1 - i/(Omega-1)), Nf)`, whose exponent never sees
#' `n_segments`, so with `N0 = 3000` and `n_segments = 200` the count of
#' effective iterations is `ceiling((Omega-1)(1 - log_N0 Nf)) + 1`: 3 is
#' bit-identical to 2, and 5 performs two removal steps rather than five.
#' Belem et al. use 2 as a *speed* setting, and their reason is specific to
#' the differential IFT they optimise, which this does not use. Measured on
#' plGeoAdaptels, delineation at 2 was the worst of the values tried.
#'
#' Every label is a single 8-connected region. Do not pass these to
#' [enforce_connectivity()]: it tests 4-connectivity, which is the adaptel
#' grower's neighbourhood, and would split each superpixel into roughly
#' twenty pieces that were never disconnected.
#'
#' @section Seeds and reproducibility:
#' With `seeds = NULL` the starting seeds are drawn with R's own RNG, so
#' `set.seed()` controls them. They will not match plGeoAdaptels, which
#' draws with NumPy: `Generator.choice` cannot be reproduced outside NumPy,
#' and reimplementing an undocumented ordering detail of a third-party
#' library is what left rHRG disagreeing with `scikit-image`'s watershed.
#'
#' The sampling is not part of the algorithm -- Belem et al. call it a free
#' choice -- so `seeds` exists to take it out of the comparison. Hand both
#' implementations the same seeds and the rest agrees exactly; that is what
#' `tools/cross_validate_against_plgeoadaptels.R` does.
#'
#' @examples
#' d <- array(runif(3 * 60 * 60, 0, 255), dim = c(3, 60, 60))
#' set.seed(1)
#' out <- sicle(d, n_segments = 40, quiet = TRUE)
#' out$n_superpixels
#' @export
sicle <- function(data, mask = NULL, n_segments = 200, n_oversampling = 3000,
                  n_iterations = 2, saliency = NULL, seeds = NULL,
                  quiet = FALSE) {
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
  n_segments <- as.integer(n_segments)
  if (is.na(n_segments) || n_segments < 1L)
    stop("n_segments must be >= 1", call. = FALSE)

  size <- as.integer(rows) * as.integer(cols)
  # Row-major, band-sequential, matching the NumPy (n_layers, size) C-order
  # buffer -- see the note in adaptels().
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
  flat[is.na(flat)] <- 0

  use_sal <- !is.null(saliency)
  if (use_sal) {
    saliency <- as.matrix(saliency)
    if (!identical(dim(saliency), c(rows, cols)))
      stop("saliency must be ", rows, "x", cols, ", got ",
           paste(dim(saliency), collapse = "x"), call. = FALSE)
    sal <- as.double(as.vector(t(saliency)))
    bad <- is.na(sal) & m == 0L
    # NA here does not raise anywhere downstream: it flows into the tree
    # mean, into the relevance, and a seed of unknown relevance would
    # outrank every real one. Nodata in a saliency raster is ordinary.
    if (any(bad))
      stop("saliency is NA at ", sum(bad), " pixel(s) the mask treats as ",
           "valid. Fill them (0 means 'no object here') or extend the mask.",
           call. = FALSE)
    sal[is.na(sal)] <- 0
  } else {
    sal <- double(0)
  }

  valid <- which(m == 0L) - 1L                    # 0-based flat indices
  if (length(valid) < n_segments)
    stop("not enough valid pixels (", length(valid), ") for ", n_segments,
         " segments", call. = FALSE)

  if (is.null(seeds)) {
    n_oversampling <- as.integer(n_oversampling)
    if (n_oversampling < n_segments) {
      # SICLE only removes seeds, so the target is unreachable from there.
      # Correcting is right; doing it silently is not.
      warning("n_oversampling=", n_oversampling, " is below n_segments=",
              n_segments, ". SICLE removes seeds, it never adds them, so the ",
              "target is unreachable from there; raising it to ",
              n_segments * 10L, ".", call. = FALSE)
      n_oversampling <- n_segments * 10L
    }
    k <- min(n_oversampling, length(valid))
    seed_idx <- as.integer(sample(valid, k))
  } else {
    seeds <- as.matrix(seeds)
    if (ncol(seeds) != 2L)
      stop("seeds must be an (n, 2) matrix of (row, col) pairs, got ",
           ncol(seeds), " column(s)", call. = FALSE)
    r <- as.integer(seeds[, 1]); cc <- as.integer(seeds[, 2])
    if (any(r < 1L | r > rows | cc < 1L | cc > cols))
      stop("a seed lies outside the raster (", rows, "x", cols, ")",
           call. = FALSE)
    seed_idx <- (r - 1L) * as.integer(cols) + (cc - 1L)
    if (any(m[seed_idx + 1L] != 0L))
      stop(sum(m[seed_idx + 1L] != 0L), " seed(s) fall on nodata pixels",
           call. = FALSE)
    if (anyDuplicated(seed_idx))
      stop("seeds contains duplicate pixels", call. = FALSE)
    if (length(seed_idx) < n_segments)
      stop(length(seed_idx), " seeds cannot produce ", n_segments,
           " segments; SICLE removes seeds, it never adds them", call. = FALSE)
  }

  n_zero <- length(seed_idx)
  if (!quiet)
    message(sprintf("  Seeds: %d initial, target %d", n_zero, n_segments))

  omega <- max(as.integer(n_iterations), 2L)
  labels <- NULL

  for (iteration in seq_len(omega) - 1L) {
    res <- .Call(C_ift_fmax, flat, m, as.integer(cols), as.integer(rows),
                 seed_idx)
    labels <- res[[1]]

    if (length(seed_idx) <= n_segments) break

    rel <- .Call(C_seed_relevance, flat, m, as.integer(cols), as.integer(rows),
                 labels, length(seed_idx), sal, use_sal)

    t <- if (omega > 1L) (iteration + 1) / (omega - 1) else 1
    m_keep <- max(as.integer(n_zero^(1 - t)), n_segments)
    m_keep <- min(m_keep, length(seed_idx))

    # Rank NaN last, explicitly. Sorting descending by reversing an ascending
    # sort puts NaN at the *head* in both R and NumPy, so a seed whose
    # relevance could not be computed would outrank every seed whose could.
    rank <- ifelse(is.nan(rel), -Inf, rel)

    # rev(order(.)) rather than order(., decreasing = TRUE): the Python does
    # argsort(...)[::-1], so among equal scores it takes the reverse of the
    # ascending order, and this reproduces that. It only matters if a tie
    # straddles the m_keep cut, which is rare -- on SNP_21_2020_1.tif the
    # boundary value is unique -- and NumPy's argsort is not stable anyway,
    # so neither side guarantees it. The cross-check in tools/ is an
    # equality test and would fail loudly if it ever happened.
    seed_idx <- seed_idx[rev(order(rank))[seq_len(m_keep)]]

    if (!quiet)
      message(sprintf("  Iteration %d/%d: %d seeds remaining",
                      iteration + 1L, omega, length(seed_idx)))
  }

  list(labels = matrix(labels, nrow = rows, ncol = cols, byrow = TRUE),
       n_superpixels = length(seed_idx))
}
