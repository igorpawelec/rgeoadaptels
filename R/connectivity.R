#' Split adaptels that are not single connected regions
#'
#' @param labels Integer matrix of adaptel ids, as returned by [adaptels()].
#'   Negative values are nodata and pass through untouched.
#' @param min_size Fragments of at most this many pixels are absorbed into an
#'   adjacent adaptel instead of becoming their own. 0 keeps every fragment.
#'
#' @return A list with `labels` and `n_adaptels`.
#'
#' @details
#' Adaptels compete for pixels: a later one takes a pixel from an earlier one
#' whenever it arrives with a smaller accumulated distance. That competition
#' is what gives the method its boundary adherence, and it can also cut an
#' earlier adaptel in two, leaving one id spread over separate patches. On a
#' 400x400 three-band raster at the default threshold, about 10 per cent of
#' adaptels come out in more than one piece.
#'
#' Whether it matters depends on what comes next. It is harmless when the
#' labels are only a lookup. It is not harmless for object-based analysis:
#' zonal statistics over a split adaptel average two spatially separate
#' patches into one "object".
#'
#' This is not applied automatically, because it changes the adaptel count.
#'
#' Uses 4-connectivity, matching the grower's default neighbourhood. If you
#' grew with `queen_topology = TRUE`, an adaptel that is 8-connected but not
#' 4-connected will be split here -- that is a real mismatch of conventions,
#' not a defect, and worth knowing before you read the count.
#'
#' @examples
#' d <- array(runif(3 * 40 * 40, 0, 255), dim = c(3, 40, 40))
#' seg <- adaptels(d, threshold = 60)
#' split <- enforce_connectivity(seg$labels)
#' c(before = seg$n_adaptels, after = split$n_adaptels)
#' @export
enforce_connectivity <- function(labels, min_size = 0) {
  if (!is.numeric(labels))
    stop("labels must be numeric", call. = FALSE)
  # Checked before as.matrix, not after: as.matrix() on a 3-D array does not
  # fail, it silently reshapes to n x 1, so a validation placed downstream of
  # it passes on exactly the input it was meant to reject.
  nd <- length(dim(labels))
  if (nd != 2L)
    stop("labels must be 2-D, got ", if (nd == 0L) "a vector" else
         paste0(nd, " dimensions"), call. = FALSE)
  labels <- as.matrix(labels)
  min_size <- as.integer(min_size)
  if (is.na(min_size) || min_size < 0L)
    stop("min_size must be >= 0", call. = FALSE)

  rows <- nrow(labels); cols <- ncol(labels)
  # t(): the C walks the buffer row-major, as the Python does.
  flat <- as.integer(as.vector(t(labels)))

  res <- .Call(C_enforce_connectivity, flat, as.integer(rows),
               as.integer(cols), min_size)

  list(labels = matrix(res[[1]], nrow = rows, ncol = cols, byrow = TRUE),
       n_adaptels = res[[2]])
}
