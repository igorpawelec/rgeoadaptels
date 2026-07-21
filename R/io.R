#' Read raster bands for segmentation
#'
#' @param x A `SpatRaster`, or a path, or a character vector of paths whose
#'   first band each becomes one layer.
#'
#' @return A list with `data`, a `(bands, rows, cols)` array; `mask`, a
#'   `(rows, cols)` integer matrix with 1 for nodata; and `template`, the
#'   `SpatRaster` the geometry came from.
#'
#' @details
#' Requires \pkg{terra}, which is only a Suggests: [adaptels()] and
#' [enforce_connectivity()] work on plain arrays and have no spatial
#' dependency at all.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   r <- terra::rast(nrows = 8, ncols = 6, nlyrs = 3,
#'                    vals = runif(144, 0, 255))
#'   b <- read_bands(r)
#'   dim(b$data)
#' }
#' }
#' @export
read_bands <- function(x) {
  .need_terra()
  r <- if (inherits(x, "SpatRaster")) x else if (length(x) > 1L)
    terra::rast(lapply(x, function(p) terra::rast(p)[[1]])) else terra::rast(x)

  nb <- terra::nlyr(r)
  first <- terra::as.matrix(r[[1]], wide = TRUE)
  rows <- nrow(first); cols <- ncol(first)
  out <- array(0, dim = c(nb, rows, cols))
  na <- matrix(FALSE, rows, cols)
  for (l in seq_len(nb)) {
    m <- terra::as.matrix(r[[l]], wide = TRUE)
    na <- na | is.na(m)
    m[is.na(m)] <- 0
    out[l, , ] <- m
  }
  list(data = out, mask = matrix(as.integer(na), rows, cols), template = r)
}


#' Segment a raster and write the labels
#'
#' @param input A `SpatRaster`, a path, or several paths taken as bands.
#' @param output Path for the label raster. `NULL` writes nothing.
#' @param ... Passed to [adaptels()].
#' @param connectivity Apply [enforce_connectivity()] to the result before
#'   writing. Off by default, because it changes the adaptel count.
#' @param min_size Passed to [enforce_connectivity()] when `connectivity` is
#'   `TRUE`.
#'
#' @return Invisibly, the list [adaptels()] returns.
#'
#' @details
#' Nodata is written as `-9999`, matching plGeoAdaptels, and set as the
#' raster's NA value so that terra and GDAL both read it back as missing
#' rather than as an adaptel numbered minus nine thousand.
#'
#' Requires \pkg{terra}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   r <- terra::rast(nrows = 20, ncols = 24, nlyrs = 3,
#'                    vals = runif(1440, 0, 255))
#'   seg <- adaptels_raster(r, output = NULL, threshold = 60)
#'   seg$n_adaptels
#' }
#' }
#' @export
adaptels_raster <- function(input, output = NULL, ..., connectivity = FALSE,
                            min_size = 0) {
  .need_terra()
  b <- read_bands(input)
  seg <- adaptels(b$data, mask = b$mask, ...)
  if (connectivity) seg <- enforce_connectivity(seg$labels, min_size = min_size)

  if (!is.null(output)) {
    r <- terra::rast(b$template[[1]])
    # t(): terra fills a layer row-major and the labels matrix is row-major
    # already, so the transpose puts them back in step. Without it the
    # output is the transpose of the input geometry, which looks plausible
    # on a square raster and is wrong on every other one.
    terra::values(r) <- as.vector(t(seg$labels))
    names(r) <- "adaptel"
    terra::NAflag(r) <- -9999
    terra::writeRaster(r, output, overwrite = TRUE, datatype = "INT4S")
  }
  invisible(seg)
}

.need_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("this function needs the 'terra' package; the segmentation itself ",
         "does not", call. = FALSE)
}
