/* Binary min-heap, shared by the adaptel grower and the SICLE forest.
 *
 * 1-based, as in the Python this is ported from, so the sift arithmetic is
 * the same expression rather than a translated one. Slot 0 is unused.
 *
 * Grows rather than caps. plGeoAdaptels capped both of its heaps and wrote
 * the pixel's label and cost *before* checking capacity, so an overflow
 * marked a pixel conquered that never propagated. That fired in SICLE --
 * 83 lost pixels on a 1200x1200 raster -- and was one insert away from
 * firing in the grower.
 *
 * Memory comes from R_alloc, which cannot be realloc'ed, so growth copies
 * into a fresh block. Everything is released when .Call returns.
 *
 * Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
 */
#ifndef RGA_HEAP_H
#define RGA_HEAP_H

#include <R.h>
#include <Rinternals.h>

typedef struct {
    double *dist;
    int *x;
    int *y;
    R_xlen_t *idx;
    R_xlen_t n;      /* elements held; slot 0 unused */
    R_xlen_t cap;    /* slots allocated, excluding slot 0 */
} Heap;

static void heap_init(Heap *h, R_xlen_t cap) {
    if (cap < 16) cap = 16;
    h->cap = cap;
    h->n = 0;
    h->dist = (double *) R_alloc(cap + 1, sizeof(double));
    h->x = (int *) R_alloc(cap + 1, sizeof(int));
    h->y = (int *) R_alloc(cap + 1, sizeof(int));
    h->idx = (R_xlen_t *) R_alloc(cap + 1, sizeof(R_xlen_t));
}

static void heap_grow(Heap *h) {
    R_xlen_t cap = h->cap * 2;
    double *d = (double *) R_alloc(cap + 1, sizeof(double));
    int *x = (int *) R_alloc(cap + 1, sizeof(int));
    int *y = (int *) R_alloc(cap + 1, sizeof(int));
    R_xlen_t *ix = (R_xlen_t *) R_alloc(cap + 1, sizeof(R_xlen_t));
    for (R_xlen_t i = 0; i <= h->n; i++) {
        d[i] = h->dist[i]; x[i] = h->x[i]; y[i] = h->y[i]; ix[i] = h->idx[i];
    }
    h->dist = d; h->x = x; h->y = y; h->idx = ix; h->cap = cap;
}

static void heap_insert(Heap *h, double dist, int x, int y, R_xlen_t idx) {
    if (h->n + 2 > h->cap) heap_grow(h);
    h->n++;
    R_xlen_t pos = h->n;
    h->dist[pos] = dist; h->x[pos] = x; h->y[pos] = y; h->idx[pos] = idx;
    while (pos > 1) {
        R_xlen_t parent = pos / 2;
        if (h->dist[pos] < h->dist[parent]) {
            double td = h->dist[pos]; h->dist[pos] = h->dist[parent]; h->dist[parent] = td;
            int tx = h->x[pos]; h->x[pos] = h->x[parent]; h->x[parent] = tx;
            int ty = h->y[pos]; h->y[pos] = h->y[parent]; h->y[parent] = ty;
            R_xlen_t ti = h->idx[pos]; h->idx[pos] = h->idx[parent]; h->idx[parent] = ti;
            pos = parent;
        } else break;
    }
}

static void heap_extract(Heap *h, double *dist, int *x, int *y, R_xlen_t *idx) {
    *dist = h->dist[1]; *x = h->x[1]; *y = h->y[1]; *idx = h->idx[1];
    h->dist[1] = h->dist[h->n];
    h->x[1] = h->x[h->n];
    h->y[1] = h->y[h->n];
    h->idx[1] = h->idx[h->n];
    h->n--;
    R_xlen_t i = 1;
    for (;;) {
        R_xlen_t smallest = i, left = i * 2, right = i * 2 + 1;
        if (left <= h->n && h->dist[left] < h->dist[smallest]) smallest = left;
        if (right <= h->n && h->dist[right] < h->dist[smallest]) smallest = right;
        if (smallest != i) {
            double td = h->dist[i]; h->dist[i] = h->dist[smallest]; h->dist[smallest] = td;
            int tx = h->x[i]; h->x[i] = h->x[smallest]; h->x[smallest] = tx;
            int ty = h->y[i]; h->y[i] = h->y[smallest]; h->y[smallest] = ty;
            R_xlen_t ti = h->idx[i]; h->idx[i] = h->idx[smallest]; h->idx[smallest] = ti;
            i = smallest;
        } else break;
    }
}

#endif /* RGA_HEAP_H */
