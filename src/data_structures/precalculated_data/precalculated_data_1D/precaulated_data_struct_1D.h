#ifndef PRECALCULATED_DATA_STRUCT_1D_H
#define PRECALCULATED_DATA_STRUCT_1D_H

#include "global_float_precision.h"

struct PrecalculatedVariables2D
{
    FLOAT_P *one_over_rho0;
    FLOAT_P *grad_T0;
    FLOAT_P *grad_g;
    FLOAT_P *grad_rho0;
    FLOAT_P *eta_over_four_pi_rho0_T0;

    FLOAT_P *VIS_COEFF_over_rho0;
    FLOAT_P *VIS_COEFF_over_T0_rho0;
    FLOAT_P *THERM_COEFF_over_T0_rho0;

    FLOAT_P two_VIS_COEFF;
    FLOAT_P one_over_2dz;
    FLOAT_P one_over_dzdz;
    FLOAT_P one_over_dz;
};

#endif // PRECALCULATED_DATA_STRUCT_1D_H