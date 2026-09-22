# The README figures, generated from test_data so they can be remade.
#
#     Rscript www/figures.R           # from the repository root
#
# Writes three figures and prints the numbers the README captions quote:
#   www/pipeline.png    the orthophoto with hand-placed dead-tree points, the
#                       adaptels at the default threshold, and the crowns
#                       grow_seeds grows from the points;
#   www/threshold.png   a 30 x 30 m window at threshold 30, 60 and 120;
#   www/max_cost.png    grow_seeds at max_cost 8, 15 and 25.
# Same scene, parameters and palette as pygeoadaptels' www/figures.py, so the
# two READMEs show the same thing from either language -- and the numbers
# printed here are expected to match the Python ones exactly, because the two
# packages are bit-identical. Uses the package in this checkout, not an
# installed copy. The orthophoto is contrast-stretched for display only.
suppressPackageStartupMessages({ library(terra) })
pkgload::load_all(".", quiet = TRUE)

rgb_file <- "test_data/SNP_21_2020_1.tif"
lab_file <- "test_data/SNP_21_2020_1_lab.tif"
pts_file <- "test_data/dead_trees_test.shp"
orange <- "#eb6834"; ink <- "#1f1f1f"
recipe <- list(band_weights = c(0.5, 2.5, 1), max_radius = 20, fill_holes = TRUE)   # docs/grow_seeds_guide.md
window <- c(row0 = 120, col0 = 90, n = 120)   # 0-based like the Python script; a 30 x 30 m window with dead trees and shadow

b  <- read_bands(rgb_file)          # data (bands, rows, cols), mask 1 = nodata, template
lb <- read_bands(lab_file)
tmpl <- b$template; res_m <- terra::res(tmpl)[1]; e <- terra::ext(tmpl)
rows <- dim(b$data)[2]; cols <- dim(b$data)[3]
valid <- b$mask == 0

# points -> 1-based (row, col) pixel pairs, the way grow_seeds_raster does it
xy <- terra::crds(terra::vect(pts_file))
seeds <- cbind(row = floor((e$ymax - xy[, 2]) / res_m) + 1, col = floor((xy[, 1] - e$xmin) / res_m) + 1)
storage.mode(seeds) <- "integer"

# ---- the runs --------------------------------------------------------------------
adapt <- lapply(c(30, 60, 120), function(t) adaptels(b$data, mask = b$mask, threshold = t))
names(adapt) <- c("30", "60", "120")
grown <- lapply(c(8, 15, 25), function(cst) do.call(grow_seeds, c(list(lb$data, seeds, mask = lb$mask, max_cost = cst, quiet = TRUE), recipe)))
names(grown) <- c("8", "15", "25")
n_px <- sum(valid)
r0 <- window[["row0"]]; c0 <- window[["col0"]]; n <- window[["n"]]
win_rows <- (r0 + 1):(r0 + n); win_cols <- (c0 + 1):(c0 + n)
cat(sprintf("scene %s: %d x %d px at %g m, %s px inside the plot, %d dead-tree points; window rows %d-%d, cols %d-%d = %g m\n",
            basename(rgb_file), cols, rows, res_m, format(n_px, big.mark = ","), nrow(seeds), r0, r0 + n, c0, c0 + n, n * res_m))
for (t in names(adapt)) {
  lab <- adapt[[t]]$labels; ok <- lab >= 0
  sizes <- tabulate(lab[ok] + 1L)
  w <- lab[win_rows, win_cols]
  cat(sprintf("threshold %5s: %s adaptels, median %.1f m2, largest %.0f m2; %d in the window\n", t,
              format(adapt[[t]]$n_adaptels, big.mark = ","), median(sizes) * res_m^2, max(sizes) * res_m^2,
              length(unique(w[w >= 0]))))
}
for (cst in names(grown)) {
  g <- grown[[cst]]$labels; ok <- g >= 0
  sizes <- tabulate(g[ok] + 1L, nbins = nrow(seeds))
  cat(sprintf("max_cost %4s: %d of %d points grew, median crown %.1f m2, largest %.0f m2, %.1f%% of the plot\n", cst,
              sum(sizes > 0), nrow(seeds), median(sizes[sizes > 0]) * res_m^2, max(sizes) * res_m^2, 100 * sum(ok) / n_px))
}

# ---- drawing ---------------------------------------------------------------------
stretched <- function() {                    # per-band percentile stretch of the valid pixels, display only
  out <- array(255, dim = dim(b$data))
  for (k in 1:3) {
    v <- b$data[k, , ][valid]; q <- quantile(v, c(0.01, 0.995))
    s <- pmin(pmax((b$data[k, , ] - q[1]) / max(q[2] - q[1], 1) * 255, 0), 255)
    s[!valid] <- 255
    out[k, , ] <- s
  }
  r <- terra::rast(tmpl, nlyrs = 3)
  terra::values(r) <- cbind(as.numeric(t(out[1, , ])), as.numeric(t(out[2, , ])), as.numeric(t(out[3, , ])))
  r
}
img <- stretched()
win_ext <- terra::ext(e$xmin + c0 * res_m, e$xmin + (c0 + n) * res_m, e$ymax - (r0 + n) * res_m, e$ymax - r0 * res_m)

layer <- function(m) { r <- terra::rast(tmpl, nlyrs = 1); terra::values(r) <- as.numeric(t(m)); r }
edges_of <- function(lab) {                  # one pixel per border, not one on each side
  ok <- lab >= 0
  edge <- matrix(FALSE, nrow(lab), ncol(lab))
  edge[, -ncol(lab)] <- lab[, -ncol(lab)] != lab[, -1]
  edge[-nrow(lab), ] <- edge[-nrow(lab), ] | (lab[-nrow(lab), ] != lab[-1, ])
  edge & ok
}
panel <- function(label, ext = e) {
  terra::plotRGB(terra::crop(img, ext), stretch = NULL, axes = FALSE, mar = c(3.2, 0.4, 0.4, 0.4), maxcell = Inf)
  rect(ext$xmin, ext$ymin, ext$xmax, ext$ymax, border = "#d9d9d9")
  h <- ext$ymax - ext$ymin
  text((ext$xmin + ext$xmax) / 2, ext$ymin - 0.05 * h, label, adj = c(0.5, 1), cex = 0.86, col = ink, xpd = NA)
}
overlay <- function(lab, fill = FALSE, alpha = 0.32, edge_alpha = 1, ext = e) {
  if (fill) {
    f <- layer(ifelse(lab >= 0, 1, NA))
    terra::plot(terra::crop(f, ext), col = adjustcolor(orange, alpha), legend = FALSE, axes = FALSE, add = TRUE, maxcell = Inf)
  }
  ed <- layer(ifelse(edges_of(lab), 1, NA))
  terra::plot(terra::crop(ed, ext), col = adjustcolor(orange, edge_alpha), legend = FALSE, axes = FALSE, add = TRUE, maxcell = Inf)
}
draw_points <- function(cex = 0.95) points(xy, pch = 21, bg = ink, col = "white", cex = cex, lwd = 0.9)

# ---- figure 1: the pipeline ----------------------------------------------------
png("www/pipeline.png", width = 11.4, height = 4.1, units = "in", res = 160, bg = "white")
layout(matrix(1:3, 1))
panel(sprintf("orthophoto, %g m, %d dead-tree points", res_m, nrow(seeds))); draw_points()
panel(sprintf("%s adaptels at threshold 60", format(adapt[["60"]]$n_adaptels, big.mark = ","))); overlay(adapt[["60"]]$labels, edge_alpha = 0.5)
panel("crowns grown from the points, max_cost 15"); overlay(grown[["15"]]$labels, fill = TRUE); draw_points(cex = 0.6)
invisible(dev.off())

# ---- figure 2: the adaptel threshold, in a window ----------------------------------
png("www/threshold.png", width = 11.4, height = 4.1, units = "in", res = 160, bg = "white")
layout(matrix(1:3, 1))
for (t in names(adapt)) {
  lab <- adapt[[t]]$labels; w <- lab[win_rows, win_cols]
  panel(sprintf("threshold = %s\n%d adaptels in this 30 m window", t, length(unique(w[w >= 0]))), ext = win_ext)
  overlay(lab, ext = win_ext)
}
invisible(dev.off())

# ---- figure 3: the grow_seeds tolerance --------------------------------------------
png("www/max_cost.png", width = 11.4, height = 4.1, units = "in", res = 160, bg = "white")
layout(matrix(1:3, 1))
for (cst in names(grown)) {
  g <- grown[[cst]]$labels
  panel(sprintf("max_cost = %s\n%.1f%% of the plot in crowns", cst, 100 * sum(g >= 0) / n_px))
  overlay(g, fill = TRUE); draw_points(cex = 0.6)
}
invisible(dev.off())
cat("written: www/pipeline.png, www/threshold.png, www/max_cost.png\n")
