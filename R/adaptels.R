DISTANCES <- c(minkowski = 0L, cosine = 1L, angular = 2L)

# Both entry points check the same things, and they check them here so they
# cannot drift apart. plGeoAdaptels shipped with adaptels_from_array()
# quietly falling back to minkowski on an unrecognised distance while
# create_adaptels() raised on the same input.
.validate <- function(threshold, distance, minkowski_p) {
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold))
    stop("threshold must be a single number", call. = FALSE)
  if (threshold <= 0)
    stop("threshold must be greater than 0", call. = FALSE)
  if (!is.character(distance) || length(distance) != 1L ||
      !distance %in% names(DISTANCES))
    stop("unknown distance '", distance, "'. Use one of: ",
         paste(names(DISTANCES), collapse = ", "), call. = FALSE)
  if (!is.numeric(minkowski_p) || length(minkowski_p) != 1L || minkowski_p <= 0)
    stop("minkowski_p must be a positive number", call. = FALSE)

  # The metrics live on different scales, so one threshold does not carry
  # across them. cosine and angular are bounded by 1 by construction; a
  # threshold above that means "merge everything" and returns a single
  # adaptel, which reads as a broken algorithm rather than a misread
  # parameter. minkowski grows with the data range and has no bound to check.
  if (distance %in% c("cosine", "angular") && threshold > 1)
    stop("threshold=", threshold, " is outside the range of the '", distance,
         "' metric, which produces distances in [0, 1]. Anything above 1 ",
         "merges the whole raster into one adaptel. Typical values are ",
         "0.005-0.2; try 0.03 to start. (The default of 60 is scaled for ",
         "'minkowski', whose distances follow the data range.)",
         call. = FALSE)
  DISTANCES[[distance]]
}


#' Create adaptels from raster bands
#'
#' Scale-adaptive superpixels: regions grow from a seed until their internal
#' distance exceeds `threshold`, and the pixels beyond it become seeds in
#' turn. Unlike SLIC and its relatives there is no target count -- the number
#' of adaptels follows the scene.
#'
#' @param data Numeric matrix `(rows, cols)` for one band, or array
#'   `(bands, rows, cols)` for several.
#' @param mask Optional matrix `(rows, cols)`. Non-zero marks nodata. When
#'   `NULL`, `NA` in any band is treated as nodata.
#' @param threshold Growth threshold. **Per metric**, not universal: 60 suits
#'   `minkowski` on raw data, while `cosine` and `angular` are bounded by 1
#'   and want something around 0.03.
#' @param distance One of `"minkowski"`, `"cosine"`, `"angular"`.
#' @param minkowski_p Exponent for `"minkowski"`. 2 is Euclidean.
#' @param queen_topology `TRUE` for 8-connectivity, `FALSE` (default) for 4.
#' @param normalize Rescale every band to [0, 1] before growing. Doing so
#'   caps the largest possible minkowski distance at `n_bands^(1/p)`, so the
#'   threshold has to be rescaled with it; passing the raw-data default is
#'   rejected rather than silently returning one adaptel.
#'
#' @return A list with `labels`, an integer matrix `(rows, cols)` holding
#'   0-based adaptel ids with `-9999` for nodata, and `n_adaptels`.
#'
#' @details
#' Adaptels compete: a later adaptel takes a pixel from an earlier one
#' whenever it reaches that pixel with a smaller accumulated distance. That
#' competition is what gives the method its boundary adherence, and it is
#' also why a label can end up in more than one piece -- see
#' [enforce_connectivity()].
#'
#' @examples
#' d <- array(runif(3 * 40 * 40, 0, 255), dim = c(3, 40, 40))
#' out <- adaptels(d, threshold = 60)
#' out$n_adaptels
#' @useDynLib rgeoadaptels, .registration = TRUE, .fixes = "C_"
#' @export
adaptels <- function(data, mask = NULL, threshold = 60,
                     distance = "minkowski", minkowski_p = 2,
                     queen_topology = FALSE, normalize = FALSE) {
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

  dist_code <- .validate(threshold, distance, minkowski_p)
  size <- as.numeric(rows) * as.numeric(cols)

  # Row-major, band-sequential: layer l, pixel (y, x) lands at
  # l*size + y*cols + x, matching the NumPy (n_layers, size) C-order buffer
  # the Python indexes. aperm puts the band axis last so that as.vector,
  # which reads column-major, walks x fastest and then y, then band.
  flat <- as.vector(aperm(arr, c(3, 2, 1)))
  flat <- as.double(flat)

  if (is.null(mask)) {
    m <- as.integer(is.na(flat[seq_len(size)]))
    for (l in seq_len(n_layers)[-1])
      m <- m | as.integer(is.na(flat[(l - 1) * size + seq_len(size)]))
    m <- as.integer(m)
  } else {
    if (!identical(dim(as.matrix(mask)), c(rows, cols)))
      stop("mask must be ", rows, "x", cols, call. = FALSE)
    m <- as.integer(as.vector(t(as.matrix(mask))) != 0)
  }
  flat[is.na(flat)] <- 0

  if (normalize) {
    ceiling_d <- n_layers^(1 / minkowski_p)
    if (distance == "minkowski" && threshold >= ceiling_d)
      stop("threshold=", threshold, " cannot be reached with normalize=TRUE. ",
           "Normalising puts every band in [0, 1], so the largest possible ",
           "minkowski distance across ", n_layers, " band(s) at p=",
           minkowski_p, " is ", format(ceiling_d, digits = 5),
           ". A threshold at or above that merges the entire raster into one ",
           "adaptel. Rescale it -- on 3-band imagery 0.1-0.5 is a working ",
           "span -- or drop normalize and keep the raw-data threshold.",
           call. = FALSE)
    for (l in seq_len(n_layers)) {
      i <- (l - 1) * size + seq_len(size)
      v <- flat[i][m == 0]
      lo <- min(v); hi <- max(v)
      if (hi > lo) flat[i] <- (flat[i] - lo) / (hi - lo) else flat[i] <- 0
    }
  }

  res <- .Call(C_create_adaptels, flat, m, as.integer(cols), as.integer(rows),
               as.double(threshold),
               as.integer(if (queen_topology) 8L else 4L),
               as.integer(dist_code), as.double(minkowski_p))

  list(labels = matrix(res[[1]], nrow = rows, ncol = cols, byrow = TRUE),
       n_adaptels = res[[2]])
}
