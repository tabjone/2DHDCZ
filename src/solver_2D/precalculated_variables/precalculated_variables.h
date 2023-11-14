#ifndef PRECALCULATED_VARIABLES_H__
#define PRECALCULATED_VARIABLES_H__

#include "global_parameters.h"
#include "global_constants.h"
#include "shared_files.h"

struct PrecalculatedVariables
{
    FLOAT_P *one_over_rho0;
    FLOAT_P *grad_T0;
    FLOAT_P *grad_g;
    FLOAT_P *grad_rho0;
    FLOAT_P *eta_over_four_pi_rho0_T0;

    FLOAT_P *VIS_COEFF_over_rho0;
    FLOAT_P *VIS_COEFF_over_T0_rho0;
    FLOAT_P *THERM_COEFF_over_T0_rho0;

    int *j_plus;
    int *j_minus;
    int *j_plus2;
    int *j_minus2;

    FLOAT_P two_VIS_COEFF;
    FLOAT_P one_over_8dydydz;
    FLOAT_P one_over_8dydzdz;
    FLOAT_P one_over_2dy;
    FLOAT_P one_over_2dz;
    FLOAT_P one_over_dydy;
    FLOAT_P one_over_dzdz;
    FLOAT_P one_over_dz;
    FLOAT_P one_over_dy;
    FLOAT_P one_over_4dydz;
};

void allocate_calculate_precalc_struct(struct PrecalculatedVariables **pv, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);
void deallocate_precalc_struct(struct PrecalculatedVariables *pv);

#endif // PRECALCULATED_VARIABLES_H__