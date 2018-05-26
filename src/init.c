#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "sglasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"sglasso_ccd_single", (DL_FUNC) &F77_SUB(sglasso_ccd_single), 24},
    {"sglasso_ccm_single", (DL_FUNC) &F77_SUB(sglasso_ccm_single), 24},
    {"sglasso_ccd_path", (DL_FUNC) &F77_SUB(sglasso_ccd_path), 26},
    {"sglasso_ccm_path", (DL_FUNC) &F77_SUB(sglasso_ccm_path), 26},
    {"gdf_fun", (DL_FUNC) &F77_SUB(gdf_fun), 8},
    {NULL, NULL, 0}
};

void attribute_visible R_init_sglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
