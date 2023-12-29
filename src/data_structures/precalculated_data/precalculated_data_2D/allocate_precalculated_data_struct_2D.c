#include <stdlib.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "precalculated_data_struct_2D.h"

void allocate_precalculate_data_struct_2D(struct PrecalculatedVariables2D **pv, int nz_full)
{
    // Allocate memory for the struct
    *pv = (struct PrecalculatedVariables2D *)malloc(sizeof(struct PrecalculatedVariables2D));
    // HANDLE MALLOC FAILURE

    // Allocate memory for the arrays
    allocate_1D_array(&(*pv)->one_over_rho0, nz_full);
    allocate_1D_array(&(*pv)->grad_g, nz_full);
    allocate_1D_array(&(*pv)->grad_rho0, nz_full);
    allocate_1D_array(&(*pv)->eta_over_four_pi_rho0_T0, nz_full);
    allocate_1D_array(&(*pv)->grad_T0, nz_full);
    allocate_1D_array(&(*pv)->VIS_COEFF_over_rho0, nz_full);
    allocate_1D_array(&(*pv)->VIS_COEFF_over_T0_rho0, nz_full);
    allocate_1D_array(&(*pv)->THERM_COEFF_over_T0_rho0, nz_full);
}