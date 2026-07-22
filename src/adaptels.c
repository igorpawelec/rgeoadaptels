/* Scale-Adaptive Superpixels (Adaptels).
 *
 * Achanta, Marquez-Neila, Fua and Susstrunk, "Scale-Adaptive Superpixels",
 * Color and Imaging Conference, 2018. Ported from the plGeoAdaptels Python
 * package, itself a reimplementation of Pawel Netzel's original C.
 *
 * Kept deliberately close to the Python, statement for statement, so the two
 * can be diffed by eye when they disagree. The details that decide whether
 * they agree at all:
 *
 *   - data arrives row-major. R passes t(matrix), whose column-major buffer
 *     is byte-identical to a row-major NumPy array, so a pixel is at
 *     y*cols + x in both. Skip that and the neighbour order transposes,
 *     which changes which adaptel claims a contested pixel;
 *   - bands are indexed layers[l * size + i], matching NumPy's (n_layers,
 *     size) C-order rather than an R matrix's column-major layout;
 *   - the heap is 1-based, as in the Python, so the sift arithmetic is the
 *     same expression and not a translated one;
 *   - both buffers grow rather than cap. plGeoAdaptels capped seeds at
 *     100000 until 0.5.0 and discarded the rest in silence, which changed
 *     27% of the partition on a 1200x1200 scene.
 *
 * Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include "heap.h"

/* ------------------------------------------------------------------ */
/* Distance between a candidate pixel and an adaptel's running colour   */
/* ------------------------------------------------------------------ */

static double calc_distance(const double *layers, int n_layers,
                            const double *cumul, R_xlen_t size, R_xlen_t index,
                            R_xlen_t sp_size_i, int dist_type, double mink_p) {
    double dist = 0.0;
    double sp = (double) sp_size_i;

    if (dist_type == 0) {                      /* minkowski */
        for (int i = 0; i < n_layers; i++) {
            double diff = fabs(cumul[i] - layers[(R_xlen_t) i * size + index] * sp);
            dist += pow(diff, mink_p);
        }
        dist = pow(dist, 1.0 / mink_p) / sp;
    } else {                                   /* cosine (1) and angular (2) */
        double A = 1.0, B = 1.0, AB = 1.0;
        for (int i = 0; i < n_layers; i++) {
            double v = layers[(R_xlen_t) i * size + index];
            double mean = cumul[i] / sp;
            A += mean * mean;
            B += v * v;
            AB += v * mean;
        }
        double val = AB / sqrt(A * B);
        if (dist_type == 1) {
            /* cosine DISTANCE, not similarity: the bare ratio is 1.0 for
             * identical spectra and falls towards 0 as they diverge, which
             * is backwards for something compared against a growth
             * threshold. plGeoAdaptels shipped it inverted until 0.3.0. */
            dist = 1.0 - val;
        } else {
            if (val > 1.0) val = 1.0;
            if (val < -1.0) val = -1.0;
            dist = acos(val) / M_PI;
        }
    }
    return dist;
}

/* ------------------------------------------------------------------ */
/* Main growth                                                          */
/* ------------------------------------------------------------------ */

SEXP C_create_adaptels(SEXP s_layers, SEXP s_mask, SEXP s_cols, SEXP s_rows,
                       SEXP s_threshold, SEXP s_connectivity,
                       SEXP s_dist_type, SEXP s_mink_p) {
    const double *layers = REAL(s_layers);
    const int *mask = INTEGER(s_mask);
    int cols = INTEGER(s_cols)[0];
    int rows = INTEGER(s_rows)[0];
    double threshold = REAL(s_threshold)[0];
    int connectivity = INTEGER(s_connectivity)[0];
    int dist_type = INTEGER(s_dist_type)[0];
    double mink_p = REAL(s_mink_p)[0];

    R_xlen_t size = (R_xlen_t) cols * (R_xlen_t) rows;
    int n_layers = (int) (XLENGTH(s_layers) / size);

    static const int DX[8] = {-1, 0, 1, 0, -1, 1, 1, -1};
    static const int DY[8] = {0, -1, 0, 1, -1, -1, 1, 1};
    R_xlen_t dIdx[8];
    dIdx[0] = -1;
    dIdx[1] = -(R_xlen_t) cols;
    dIdx[2] = 1;
    dIdx[3] = (R_xlen_t) cols;
    dIdx[4] = -1 - (R_xlen_t) cols;
    dIdx[5] = 1 - (R_xlen_t) cols;
    dIdx[6] = 1 + (R_xlen_t) cols;
    dIdx[7] = -1 + (R_xlen_t) cols;

    SEXP s_labels = PROTECT(allocVector(INTSXP, size));
    int *labels = INTEGER(s_labels);
    double *distances = (double *) R_alloc(size, sizeof(double));
    for (R_xlen_t i = 0; i < size; i++) { labels[i] = -1; distances[i] = 0.0; }

    /* Seed buffer, grown on demand. See the header note. */
    R_xlen_t seeds_cap = 4096, n_seeds = 0;
    int *seeds_x = (int *) R_alloc(seeds_cap, sizeof(int));
    int *seeds_y = (int *) R_alloc(seeds_cap, sizeof(int));
    R_xlen_t *seeds_idx = (R_xlen_t *) R_alloc(seeds_cap, sizeof(R_xlen_t));

    Heap heap;
    heap_init(&heap, 4096);

    double *cumul = (double *) R_alloc(n_layers, sizeof(double));
    int current_label = 0;
    R_xlen_t start_idx = 0;

    while (start_idx < size) {
        int found = 0;
        while (start_idx < size) {
            if (mask[start_idx] == 0 && labels[start_idx] == -1) { found = 1; break; }
            start_idx++;
        }
        if (!found) break;

        n_seeds = 1;
        seeds_x[0] = (int) (start_idx % cols);
        seeds_y[0] = (int) (start_idx / cols);
        seeds_idx[0] = start_idx;

        for (R_xlen_t si = 0; si < n_seeds; si++) {
            int s_x = seeds_x[si], s_y = seeds_y[si];
            R_xlen_t s_idx = seeds_idx[si];
            if (labels[s_idx] >= 0) continue;

            heap.n = 0;
            heap_insert(&heap, 0.0, s_x, s_y, s_idx);
            distances[s_idx] = 0.0;
            labels[s_idx] = current_label;

            R_xlen_t sp_size = 1;
            for (int l = 0; l < n_layers; l++)
                cumul[l] = layers[(R_xlen_t) l * size + s_idx];

            while (heap.n > 0) {
                double cell_dist; int cell_x, cell_y; R_xlen_t cell_idx;
                heap_extract(&heap, &cell_dist, &cell_x, &cell_y, &cell_idx);

                for (int conn = 0; conn < connectivity; conn++) {
                    int nx = cell_x + DX[conn];
                    int ny = cell_y + DY[conn];
                    R_xlen_t nidx = cell_idx + dIdx[conn];

                    if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
                    if (mask[nidx] == 1) continue;

                    if (distances[cell_idx] < threshold) {
                        if (labels[nidx] != labels[cell_idx]) {
                            double d = distances[cell_idx] +
                                calc_distance(layers, n_layers, cumul, size,
                                              nidx, sp_size, dist_type, mink_p);
                            if (d < distances[nidx] || labels[nidx] < 0) {
                                distances[nidx] = d;
                                labels[nidx] = labels[cell_idx];
                                for (int l = 0; l < n_layers; l++)
                                    cumul[l] += layers[(R_xlen_t) l * size + nidx];
                                sp_size++;
                                heap_insert(&heap, distances[nidx], nx, ny, nidx);
                            }
                        }
                    } else if (labels[nidx] < 0) {
                        if (n_seeds >= seeds_cap) {
                            R_xlen_t cap = seeds_cap * 2;
                            int *nx_ = (int *) R_alloc(cap, sizeof(int));
                            int *ny_ = (int *) R_alloc(cap, sizeof(int));
                            R_xlen_t *ni_ = (R_xlen_t *) R_alloc(cap, sizeof(R_xlen_t));
                            for (R_xlen_t k = 0; k < n_seeds; k++) {
                                nx_[k] = seeds_x[k]; ny_[k] = seeds_y[k]; ni_[k] = seeds_idx[k];
                            }
                            seeds_x = nx_; seeds_y = ny_; seeds_idx = ni_;
                            seeds_cap = cap;
                        }
                        seeds_x[n_seeds] = nx;
                        seeds_y[n_seeds] = ny;
                        seeds_idx[n_seeds] = nidx;
                        n_seeds++;
                    }
                }
            }
            current_label++;
        }
        start_idx++;
    }

    /* Renumber consecutively in raster order; nodata becomes -9999. */
    int max_lab = current_label > 0 ? current_label : 1;
    int *lut = (int *) R_alloc(max_lab, sizeof(int));
    for (int i = 0; i < max_lab; i++) lut[i] = -1;

    for (R_xlen_t i = 0; i < size; i++) if (mask[i] == 1) labels[i] = -9999;

    int new_id = 0;
    for (R_xlen_t i = 0; i < size; i++)
        if (mask[i] == 0 && labels[i] >= 0 && lut[labels[i]] == -1)
            lut[labels[i]] = new_id++;
    for (R_xlen_t i = 0; i < size; i++)
        if (mask[i] == 0 && labels[i] >= 0) labels[i] = lut[labels[i]];

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, s_labels);
    SET_VECTOR_ELT(out, 1, ScalarInteger(new_id));
    UNPROTECT(2);
    return out;
}
