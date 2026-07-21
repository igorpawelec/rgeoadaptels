/* Split adaptels that are not single connected regions.
 *
 * Ported from plGeoAdaptels, which builds this on scipy.ndimage. Three
 * details decide whether the two produce the same ids, and none of them is
 * visible from the docstring:
 *
 *   - components are numbered in raster-scan order of first encounter,
 *     which is what ndimage.label does;
 *   - adaptels are visited in ascending id order, and new ids are handed
 *     out in that order, then by component number;
 *   - a sliver is absorbed into whichever labelled neighbour it touches
 *     most, and np.bincount(...).argmax() breaks a tie towards the
 *     *smallest* id. The neighbours form a set, not a multiset: the Python
 *     dilates to a boolean ring, so a pixel adjacent to three fragment
 *     pixels still counts once.
 *
 * 4-connectivity throughout, matching the neighbourhood the grower uses by
 * default.
 *
 * Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

static const int NDX[4] = {-1, 1, 0, 0};
static const int NDY[4] = {0, 0, -1, 1};

SEXP C_enforce_connectivity(SEXP s_labels, SEXP s_rows, SEXP s_cols,
                            SEXP s_min_size) {
    const int *labels = INTEGER(s_labels);
    int rows = INTEGER(s_rows)[0];
    int cols = INTEGER(s_cols)[0];
    int min_size = INTEGER(s_min_size)[0];
    R_xlen_t size = (R_xlen_t) rows * (R_xlen_t) cols;

    int max_lab = -1;
    for (R_xlen_t i = 0; i < size; i++) if (labels[i] > max_lab) max_lab = labels[i];

    SEXP s_out = PROTECT(allocVector(INTSXP, size));
    int *out = INTEGER(s_out);
    for (R_xlen_t i = 0; i < size; i++) out[i] = -1;

    if (max_lab < 0) {
        SEXP res = PROTECT(allocVector(VECSXP, 2));
        SET_VECTOR_ELT(res, 0, s_out);
        SET_VECTOR_ELT(res, 1, ScalarInteger(0));
        UNPROTECT(2);
        return res;
    }

    int n_lab = max_lab + 1;
    int *r0 = (int *) R_alloc(n_lab, sizeof(int));
    int *r1 = (int *) R_alloc(n_lab, sizeof(int));
    int *c0 = (int *) R_alloc(n_lab, sizeof(int));
    int *c1 = (int *) R_alloc(n_lab, sizeof(int));
    for (int l = 0; l < n_lab; l++) { r0[l] = rows; r1[l] = -1; c0[l] = cols; c1[l] = -1; }

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            int lab = labels[(R_xlen_t) y * cols + x];
            if (lab < 0) continue;
            if (y < r0[lab]) r0[lab] = y;
            if (y > r1[lab]) r1[lab] = y;
            if (x < c0[lab]) c0[lab] = x;
            if (x > c1[lab]) c1[lab] = x;
        }
    }

    int *comp = (int *) R_alloc(size, sizeof(int));
    for (R_xlen_t i = 0; i < size; i++) comp[i] = 0;
    R_xlen_t *stack = (R_xlen_t *) R_alloc(size, sizeof(R_xlen_t));

    /* Deferred slivers, stored as one flat pixel list with offsets. */
    R_xlen_t small_cap = 1024, n_small_px = 0, n_small = 0;
    R_xlen_t *small_px = (R_xlen_t *) R_alloc(small_cap, sizeof(R_xlen_t));
    R_xlen_t off_cap = 256;
    R_xlen_t *small_off = (R_xlen_t *) R_alloc(off_cap + 1, sizeof(R_xlen_t));
    small_off[0] = 0;

    int new_id = 0;

    R_xlen_t coff_cap = 256;
    R_xlen_t *comp_off = (R_xlen_t *) R_alloc(coff_cap + 1, sizeof(R_xlen_t));

    for (int lab = 0; lab < n_lab; lab++) {
        if (r1[lab] < 0) continue;                     /* id unused */

        /* Flood every component of this label first, packed end to end on
         * the stack, because the min_size rule depends on how many there
         * are: the Python takes an early `continue` for a single-component
         * adaptel, so a small one is kept whole rather than absorbed. */
        R_xlen_t ncomp = 0, total = 0;
        comp_off[0] = 0;
        for (int y = r0[lab]; y <= r1[lab]; y++) {
            for (int x = c0[lab]; x <= c1[lab]; x++) {
                R_xlen_t seed = (R_xlen_t) y * cols + x;
                if (labels[seed] != lab || comp[seed] != 0) continue;

                R_xlen_t first = total;
                stack[total++] = seed;
                comp[seed] = 1;
                while (first < total) {
                    R_xlen_t cur = stack[first++];
                    int cy = (int) (cur / cols), cx = (int) (cur % cols);
                    for (int k = 0; k < 4; k++) {
                        int nx = cx + NDX[k], ny = cy + NDY[k];
                        if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
                        R_xlen_t nidx = (R_xlen_t) ny * cols + nx;
                        if (labels[nidx] != lab || comp[nidx] != 0) continue;
                        comp[nidx] = 1;
                        stack[total++] = nidx;
                    }
                }
                if (ncomp + 1 > coff_cap) {
                    R_xlen_t cap = coff_cap * 2;
                    R_xlen_t *p = (R_xlen_t *) R_alloc(cap + 1, sizeof(R_xlen_t));
                    for (R_xlen_t i = 0; i <= ncomp; i++) p[i] = comp_off[i];
                    comp_off = p; coff_cap = cap;
                }
                comp_off[++ncomp] = total;
            }
        }

        for (R_xlen_t c = 0; c < ncomp; c++) {
            R_xlen_t from = comp_off[c], to = comp_off[c + 1], len = to - from;
            int defer = (ncomp > 1) && min_size > 0 && len <= (R_xlen_t) min_size;
            if (defer) {
                while (n_small_px + len > small_cap) {
                    R_xlen_t cap = small_cap * 2;
                    R_xlen_t *p = (R_xlen_t *) R_alloc(cap, sizeof(R_xlen_t));
                    for (R_xlen_t i = 0; i < n_small_px; i++) p[i] = small_px[i];
                    small_px = p; small_cap = cap;
                }
                if (n_small + 1 > off_cap) {
                    R_xlen_t cap = off_cap * 2;
                    R_xlen_t *p = (R_xlen_t *) R_alloc(cap + 1, sizeof(R_xlen_t));
                    for (R_xlen_t i = 0; i <= n_small; i++) p[i] = small_off[i];
                    small_off = p; off_cap = cap;
                }
                for (R_xlen_t i = from; i < to; i++) small_px[n_small_px++] = stack[i];
                small_off[++n_small] = n_small_px;
            } else {
                for (R_xlen_t i = from; i < to; i++) out[stack[i]] = new_id;
                new_id++;
            }
        }
    }

    /* Absorb slivers. Done after the main pass so they can attach to ids
     * created during it, exactly as the Python does. */
    int *ring_mark = (int *) R_alloc(size, sizeof(int));
    for (R_xlen_t i = 0; i < size; i++) ring_mark[i] = 0;
    int *cand = (int *) R_alloc(size, sizeof(int));

    for (R_xlen_t f = 0; f < n_small; f++) {
        R_xlen_t from = small_off[f], to = small_off[f + 1];
        R_xlen_t n_cand = 0;
        for (R_xlen_t i = from; i < to; i++) {
            R_xlen_t p = small_px[i];
            int py = (int) (p / cols), px = (int) (p % cols);
            for (int k = 0; k < 4; k++) {
                int nx = px + NDX[k], ny = py + NDY[k];
                if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
                R_xlen_t nidx = (R_xlen_t) ny * cols + nx;
                if (ring_mark[nidx] == 1) continue;      /* set, not multiset */
                int in_frag = 0;
                for (R_xlen_t j = from; j < to; j++)
                    if (small_px[j] == nidx) { in_frag = 1; break; }
                if (in_frag) continue;
                ring_mark[nidx] = 1;
                if (out[nidx] >= 0) cand[n_cand++] = out[nidx];
            }
        }
        for (R_xlen_t i = from; i < to; i++) {           /* clear the marks */
            R_xlen_t p = small_px[i];
            int py = (int) (p / cols), px = (int) (p % cols);
            for (int k = 0; k < 4; k++) {
                int nx = px + NDX[k], ny = py + NDY[k];
                if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) continue;
                ring_mark[(R_xlen_t) ny * cols + nx] = 0;
            }
        }

        int winner;
        if (n_cand > 0) {
            /* mode, ties to the smallest id, matching bincount().argmax() */
            int best = -1; R_xlen_t best_n = 0;
            for (R_xlen_t i = 0; i < n_cand; i++) {
                R_xlen_t c = 0;
                for (R_xlen_t j = 0; j < n_cand; j++) if (cand[j] == cand[i]) c++;
                if (c > best_n || (c == best_n && cand[i] < best)) {
                    best_n = c; best = cand[i];
                }
            }
            winner = best;
        } else {
            winner = new_id++;                  /* nothing adjacent: keep it */
        }
        for (R_xlen_t i = from; i < to; i++) out[small_px[i]] = winner;
    }

    /* Nodata keeps whatever negative marker it arrived with — -9999 from
     * the grower, not the -1 this buffer was initialised to. Matches
     * `out[~valid] = labels[~valid]`. */
    for (R_xlen_t i = 0; i < size; i++) if (labels[i] < 0) out[i] = labels[i];

    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, s_out);
    SET_VECTOR_ELT(res, 1, ScalarInteger(new_id));
    UNPROTECT(2);
    return res;
}
