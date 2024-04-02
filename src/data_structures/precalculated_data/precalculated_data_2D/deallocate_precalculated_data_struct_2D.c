#include <stdlib.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "precalculated_data_struct_2D.h"

void deallocate_precalculated_data_struct_2D(struct PrecalculatedVariables2D **pv)
{
    free((*pv)->one_over_rho0);
    free((*pv)->grad_g);
    free((*pv)->grad_rho0);
    free((*pv)->eta_over_four_pi_rho0_T0);
    free((*pv)->grad_T0);
    free((*pv)->VIS_COEFF_over_rho0);
    free((*pv)->VIS_COEFF_over_T0_rho0);
    free((*pv)->THERM_COEFF_over_T0_rho0);
    free((*pv)->damping_factor);
    free((*pv)->vz_horizontal_average);
    free(*pv);
}