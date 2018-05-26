#include <R_ext/RS.h>

void
F77_SUB(sglasso_ccd_single)(int *nSv, double *Sv, int *nTv, int *Tv_pkg,
                            int *Tv_rw, int *nv, int *Tv_ptr, int *nTe,
                            int *Te, int *nTe_ptr, int *Te_ptr, int *ne,
                            int *Te_scn, int *Te_ptr_scn, int *nstep,
                            double *trnc, double *tol, double *rho,
                            double *grd, double *th, double *w, int *pnl_flg,
                            int *n, int *conv);

void
F77_SUB(sglasso_ccm_single)(int *nSv, double *Sv, int *nTv, int *Tv_pkg,
                            int *Tv_rw, int *nv, int *Tv_ptr, int *nTe,
                            int *Te, int *nTe_ptr, int *Te_ptr, int *ne,
                            int *Te_scn, int *Te_ptr_scn, int *nstep,
                            double *trnc, double *tol, double *rho,
                            double *grd, double *th, double *w, int *pnl_flg,
                            int *n, int *conv);

void
F77_SUB(sglasso_ccd_path)(int *nSv, double *Sv, int *nTv, int *Tv_pkg,
                          int *Tv_rw, int *nv, int *Tv_ptr, int *nTe,
                          int *Te, int *nTe_ptr, int *Te_ptr, int *ne,
                          int *Te_scn, int *Te_ptr_scn, int *nstep,
                          double *trnc, double *tol, double *rho,
                          int *nrho, double *min_rho, double *grd,
                          double *th, double *w, int *pnl_flg, int *n,
                          int *conv);

void
F77_SUB(sglasso_ccm_path)(int *nSv, double *Sv, int *nTv, int *Tv_pkg,
                          int *Tv_rw, int *nv, int *Tv_ptr, int *nTe,
                          int *Te, int *nTe_ptr, int *Te_ptr, int *ne,
                          int *Te_scn, int *Te_ptr_scn, int *nstep,
                          double *trnc, double *tol, double *rho,
                          int *nrho, double *min_rho, double *grd,
                          double *th, double *w, int *pnl_flg, int *n,
                          int *conv);

void
F77_SUB(gdf_fun)(int *n, int *p, double *X, double *S, int *nrho, double *Kh,
                 double *gdf, int *info);
