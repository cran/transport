#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    _transport_SolveHierarchicalTransport

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void aha_compute_transport(void *, void *, void *, void *, void *, void *);
extern void aha_dphi(void *, void *, void *, void *, void *, void *, void *, void *);
extern void aha_free();
extern void aha_get_transport(void *, void *, void *, void *);
extern void aha_init(void *, void *, void *);
extern void aha_phi(void *, void *, void *, void *, void *, void *, void *, void *);
extern void aha_wasserstein(void *, void *, void *, void *, void *, void *);
extern void auction(void *, void *, void *, void *, void *, void *);
extern void auctionbf(void *, void *, void *, void *, void *, void *, void *);
extern void compute_power_diagram(void *, void *, void *, void *, void *, void *);
extern void decompose_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_power_diagram(void *, void *, void *);
extern void primaldual(void *, void *, void *, void *, void *, void *);
extern void revsimplex(void *, void *, void *, void *, void *, void *, void *, void *);
extern void shortsimplex(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _transport_SolveHierarchicalTransport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _transport_cplex_present();

static const R_CMethodDef CEntries[] = {
    {"aha_compute_transport", (DL_FUNC) &aha_compute_transport,  6},
    {"aha_dphi",              (DL_FUNC) &aha_dphi,               8},
    {"aha_free",              (DL_FUNC) &aha_free,               0},
    {"aha_get_transport",     (DL_FUNC) &aha_get_transport,      4},
    {"aha_init",              (DL_FUNC) &aha_init,               3},
    {"aha_phi",               (DL_FUNC) &aha_phi,                8},
    {"aha_wasserstein",       (DL_FUNC) &aha_wasserstein,        6},
    {"auction",               (DL_FUNC) &auction,                6},
    {"auctionbf",             (DL_FUNC) &auctionbf,              7},
    {"compute_power_diagram", (DL_FUNC) &compute_power_diagram,  6},
    {"decompose_c",           (DL_FUNC) &decompose_c,           10},
    {"get_power_diagram",     (DL_FUNC) &get_power_diagram,      3},
    {"primaldual",            (DL_FUNC) &primaldual,             6},
    {"revsimplex",            (DL_FUNC) &revsimplex,             8},
    {"shortsimplex",          (DL_FUNC) &shortsimplex,          10},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_transport_SolveHierarchicalTransport", (DL_FUNC) &_transport_SolveHierarchicalTransport, 13},
    {"_transport_cplex_present", (DL_FUNC) &_transport_cplex_present, 0},
    {NULL, NULL, 0}
};

void R_init_transport(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
