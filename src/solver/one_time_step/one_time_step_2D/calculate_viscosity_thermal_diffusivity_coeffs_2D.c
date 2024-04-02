#include "global_float_precision.h"
#include "global_parameters.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "spacial_derivatives_module/derivatives_2D/derivatives_2D.h"

#include <math.h>

void calculate_viscosity_thermal_diffusivity_coeffs_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j)
{
    int ny = grid_info->ny;

    // Calculating viscosity coefficient
    FLOAT_P dvy_dy = central_first_derivative_y_2D(fg->vy, i, j, ny, precalc->one_over_2dy);
    FLOAT_P dvz_dz = central_first_derivative_z_2D(fg->vz, i, j, precalc->one_over_2dz);

    FLOAT_P viscosity_coefficient = FIRST_VISCOSITY_COEFF + SECOND_VISCOSITY_COEFF*sqrt(dvy_dy*dvy_dy + dvz_dz*dvz_dz);

    precalc->two_VIS_COEFF = 2.0*viscosity_coefficient;
    precalc->VIS_COEFF_over_rho0[i] = viscosity_coefficient/bg->rho0[i];
    precalc->VIS_COEFF_over_T0_rho0[i] = viscosity_coefficient/(bg->T0[i]*bg->rho0[i]);


    // Calculating thermal diffusion coefficient
    FLOAT_P ds1_dy = central_first_derivative_y_2D(fg->s1, i, j, ny, precalc->one_over_2dy);
    FLOAT_P ds1_dz = central_first_derivative_z_2D(fg->s1, i, j, precalc->one_over_2dz);

    FLOAT_P thermal_diffusion_coefficient = FIRST_THERMAL_DIFFUSIVITY_COEFF + SECOND_THERMAL_DIFFUSIVITY_COEFF*sqrt(ds1_dy*ds1_dy + ds1_dz*ds1_dz);

    precalc->THERM_COEFF_over_T0_rho0[i] = thermal_diffusion_coefficient/(bg->T0[i]*bg->rho0[i]);
}