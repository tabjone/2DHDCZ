#include "precalculated_variables.h"
#include <stdlib.h>

void deallocate_precalc_struct(struct PrecalculatedVariables *pv)
{
    deallocate_1D_array(pv->one_over_rho0);
    deallocate_1D_array(pv->grad_g);
    deallocate_1D_array(pv->grad_rho0);
    deallocate_1D_array(pv->eta_over_four_pi_rho0_T0);
    deallocate_1D_array(pv->grad_T0);
    deallocate_1D_array(pv->VIS_COEFF_over_rho0);
    deallocate_1D_array(pv->VIS_COEFF_over_T0_rho0);
    deallocate_1D_array(pv->THERM_COEFF_over_T0_rho0);
    free(pv->j_plus);
    free(pv->j_minus);
    free(pv->j_plus2);
    free(pv->j_minus2);
    pv->j_plus = NULL;
    pv->j_minus = NULL;
    pv->j_plus2 = NULL;
    pv->j_minus2 = NULL;

    free(pv);
}