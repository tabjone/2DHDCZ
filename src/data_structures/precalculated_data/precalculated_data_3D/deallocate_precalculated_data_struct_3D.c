#include <stdlib.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "precalculated_data_struct_3D.h"

void deallocate_precalculated_data_struct_3D(struct PrecalculatedVariables3D **pv)
{
    free((*pv)->one_over_rho0);
    free((*pv)->grad_g);
    free((*pv)->grad_rho0);
    free((*pv)->eta_over_four_pi_rho0_T0);
    free((*pv)->grad_T0);
    free((*pv)->VIS_COEFF_over_rho0);
    free((*pv)->VIS_COEFF_over_T0_rho0);
    free((*pv)->THERM_COEFF_over_T0_rho0);
    free(*pv);
}