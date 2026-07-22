/* SICLE -- Superpixels through Iterative CLEarcutting.
 *
 * Belem, Barcelos, Joao, Perret, Cousty, Guimaraes and Falcao, "Novel
 * Arc-Cost Functions and Seed Relevance Estimations for Compact and
 * Accurate Superpixels", Journal of Mathematical Imaging and Vision,
 * 65:770-786, 2023. DOI 10.1007/s10851-023-01156-9
 *
 * Two kernels live here; the iteration that drives them -- the seed
 * preservation curve and the relevance ranking -- is plain R, mirroring
 * how plGeoAdaptels leaves that part to Python. It runs a handful of times
 * per segmentation, not once per pixel.
 *
 * The seeds arrive from the caller and are never drawn here. NumPy's
 * Generator.choice cannot be reproduced outside NumPy, and reimplementing
 * an undocumented ordering detail of a third-party library is what left
 * rHRG disagreeing with scikit-image's watershed on 0.25% of pixels.
 * Belem et al. treat the sampling as a free choice, so it is not part of
 * the algorithm; keeping it out of here is what lets the algorithm be
 * compared against plGeoAdaptels for equality.
 *
 * Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "heap.h"

/* 8-adjacency, in the order the Python visits it. The order decides which
 * seed claims a pixel that two reach at the same cost. */
static const int SDX[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static const int SDY[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

/* ------------------------------------------------------------------ */
/* Seed-restricted IFT, fmax path-cost, wroot arc-cost                  */
/* ------------------------------------------------------------------ */

SEXP C_ift_fmax(SEXP s_layers, SEXP s_mask, SEXP s_cols, SEXP s_rows,
                SEXP s_seeds) {
    const double *layers = REAL(s_layers);
    const int *mask = INTEGER(s_mask);
    int cols = INTEGER(s_cols)[0];
    int rows = INTEGER(s_rows)[0];
    const int *seeds = INTEGER(s_seeds);
    R_xlen_t n_seeds = XLENGTH(s_seeds);

    R_xlen_t size = (R_xlen_t) cols * (R_xlen_t) rows;
    int n_layers = (int) (XLENGTH(s_layers) / size);

    SEXP s_labels = PROTECT(allocVector(INTSXP, size));
    SEXP s_cost = PROTECT(allocVector(REALSXP, size));
    int *labels = INTEGER(s_labels);
    double *cost = REAL(s_cost);
    for (R_xlen_t i = 0; i < size; i++) { labels[i] = -1; cost[i] = 1e30; }

    /* wroot uses the seed's own features, which never change as the tree
     * grows -- that is what makes it more stable than the dynamic arc-cost
     * the paper also evaluates. Cached so the inner loop does not chase
     * the seed index through the label array. */
    double *sf = (double *) R_alloc(n_seeds * n_layers, sizeof(double));
    for (R_xlen_t si = 0; si < n_seeds; si++)
        for (int l = 0; l < n_layers; l++)
            sf[si * n_layers + l] = layers[(R_xlen_t) l * size + seeds[si]];

    Heap heap;
    heap_init(&heap, n_seeds * 2 > 4096 ? n_seeds * 2 : 4096);

    for (R_xlen_t si = 0; si < n_seeds; si++) {
        R_xlen_t idx = seeds[si];
        if (mask[idx] != 0) continue;
        cost[idx] = 0.0;
        labels[idx] = (int) si;
        heap_insert(&heap, 0.0, (int) (idx % cols), (int) (idx / cols), idx);
    }

    while (heap.n > 0) {
        double c_dist; int c_x, c_y; R_xlen_t c_idx;
        heap_extract(&heap, &c_dist, &c_x, &c_y, &c_idx);

        /* Stale entry: this pixel was reached again more cheaply after it
         * was queued. The heap keeps both; only the cheaper one counts. */
        if (c_dist > cost[c_idx]) continue;

        int seed_label = labels[c_idx];
        const double *fs = sf + (R_xlen_t) seed_label * n_layers;

        for (int k = 0; k < 8; k++) {
            int nx = c_x + SDX[k], ny = c_y + SDY[k];
            if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
            R_xlen_t nidx = (R_xlen_t) ny * cols + nx;
            if (mask[nidx] != 0) continue;

            double arc = 0.0;
            for (int l = 0; l < n_layers; l++) {
                double diff = fs[l] - layers[(R_xlen_t) l * size + nidx];
                arc += diff * diff;
            }
            arc = sqrt(arc);

            double new_cost = c_dist > arc ? c_dist : arc;   /* fmax */
            if (new_cost < cost[nidx]) {
                cost[nidx] = new_cost;
                labels[nidx] = seed_label;
                heap_insert(&heap, new_cost, nx, ny, nidx);
            }
        }
    }

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, s_labels);
    SET_VECTOR_ELT(out, 1, s_cost);
    UNPROTECT(3);
    return out;
}

/* ------------------------------------------------------------------ */
/* Seed relevance                                                       */
/* ------------------------------------------------------------------ */

SEXP C_seed_relevance(SEXP s_layers, SEXP s_mask, SEXP s_cols, SEXP s_rows,
                      SEXP s_labels, SEXP s_n_seeds, SEXP s_saliency,
                      SEXP s_use_saliency) {
    const double *layers = REAL(s_layers);
    const int *mask = INTEGER(s_mask);
    int cols = INTEGER(s_cols)[0];
    int rows = INTEGER(s_rows)[0];
    const int *labels = INTEGER(s_labels);
    R_xlen_t n_seeds = (R_xlen_t) INTEGER(s_n_seeds)[0];
    const double *sal = REAL(s_saliency);
    int use_sal = LOGICAL(s_use_saliency)[0];

    R_xlen_t size = (R_xlen_t) cols * (R_xlen_t) rows;
    int n_layers = (int) (XLENGTH(s_layers) / size);

    double *tree_sum = (double *) R_alloc(n_seeds * n_layers, sizeof(double));
    R_xlen_t *tree_count = (R_xlen_t *) R_alloc(n_seeds, sizeof(R_xlen_t));
    double *tree_sal = (double *) R_alloc(n_seeds, sizeof(double));
    for (R_xlen_t s = 0; s < n_seeds; s++) {
        tree_count[s] = 0; tree_sal[s] = 0.0;
        for (int l = 0; l < n_layers; l++) tree_sum[s * n_layers + l] = 0.0;
    }

    for (R_xlen_t i = 0; i < size; i++) {
        int lab = labels[i];
        if (lab < 0 || mask[i] != 0) continue;
        tree_count[lab]++;
        for (int l = 0; l < n_layers; l++)
            tree_sum[(R_xlen_t) lab * n_layers + l] += layers[(R_xlen_t) l * size + i];
        if (use_sal) tree_sal[lab] += sal[i];
    }

    double *tree_mean = (double *) R_alloc(n_seeds * n_layers, sizeof(double));
    double *sal_mean = (double *) R_alloc(n_seeds, sizeof(double));
    R_xlen_t total_valid = 0;
    for (R_xlen_t s = 0; s < n_seeds; s++) {
        sal_mean[s] = 0.0;
        for (int l = 0; l < n_layers; l++) tree_mean[s * n_layers + l] = 0.0;
        if (tree_count[s] > 0) {
            for (int l = 0; l < n_layers; l++)
                tree_mean[s * n_layers + l] =
                    tree_sum[s * n_layers + l] / (double) tree_count[s];
            if (use_sal) sal_mean[s] = tree_sal[s] / (double) tree_count[s];
        }
        total_valid += tree_count[s];
    }
    if (total_valid == 0) total_valid = 1;

    double *min_contrast = (double *) R_alloc(n_seeds, sizeof(double));
    double *max_sal_contrast = (double *) R_alloc(n_seeds, sizeof(double));
    for (R_xlen_t s = 0; s < n_seeds; s++) {
        min_contrast[s] = 1e30; max_sal_contrast[s] = 0.0;
    }

    /* 8-adjacency, matching the forest. Belem et al. define tree adjacency
     * over the same arc set the IFT grows on. plGeoAdaptels scanned four
     * neighbours here until 0.4.0, so a tree touching its only neighbour
     * diagonally was found to have none: its minimum contrast stayed at the
     * sentinel, collapsed to 0 below, and the seed was removed first on no
     * evidence. */
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            R_xlen_t idx = (R_xlen_t) y * cols + x;
            int lab_s = labels[idx];
            if (lab_s < 0) continue;
            const double *ms = tree_mean + (R_xlen_t) lab_s * n_layers;
            for (int k = 0; k < 8; k++) {
                int nx = x + SDX[k], ny = y + SDY[k];
                if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
                int lab_t = labels[(R_xlen_t) ny * cols + nx];
                if (lab_t < 0 || lab_t == lab_s) continue;

                const double *mt = tree_mean + (R_xlen_t) lab_t * n_layers;
                double d = 0.0;
                for (int l = 0; l < n_layers; l++) {
                    double diff = ms[l] - mt[l];
                    d += diff * diff;
                }
                d = sqrt(d);
                if (d < min_contrast[lab_s]) min_contrast[lab_s] = d;

                if (use_sal) {
                    double sd = fabs(sal_mean[lab_s] - sal_mean[lab_t]);
                    if (sd > max_sal_contrast[lab_s]) max_sal_contrast[lab_s] = sd;
                }
            }
        }
    }

    SEXP s_out = PROTECT(allocVector(REALSXP, n_seeds));
    double *rel = REAL(s_out);
    for (R_xlen_t s = 0; s < n_seeds; s++) {
        if (tree_count[s] == 0) { rel[s] = 0.0; continue; }
        double v_size = (double) tree_count[s] / (double) total_valid;
        double mc = min_contrast[s];
        if (mc > 1e29) mc = 0.0;         /* no neighbour found */
        double vsc = v_size * mc;
        rel[s] = use_sal
            ? vsc * (sal_mean[s] > max_sal_contrast[s]
                     ? sal_mean[s] : max_sal_contrast[s])
            : vsc;
    }

    UNPROTECT(1);
    return s_out;
}
