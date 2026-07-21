/* Symbol registration.
 *
 * Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_create_adaptels(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_enforce_connectivity(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"create_adaptels", (DL_FUNC) &C_create_adaptels, 8},
    {"enforce_connectivity", (DL_FUNC) &C_enforce_connectivity, 4},
    {NULL, NULL, 0}
};

void attribute_visible R_init_rgeoadaptels(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
