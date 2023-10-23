#ifndef SPACIAL_DERIVATIVES_H__
#define SPACIAL_DERIVATIVES_H__

#include "global_parameters.h"
#include "shared_files.h"
#include "../functions.h"

FLOAT_P upwind_first_derivative_y(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P dy, int ny);
FLOAT_P upwind_first_derivative_z(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P dz, int nz);

FLOAT_P central_first_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny);
FLOAT_P central_first_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz);

FLOAT_P central_second_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny);
FLOAT_P central_second_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz);

struct CommonCentralDerivatives
{
    FLOAT_P dvy_dy;
    FLOAT_P dvy_dz;
    FLOAT_P dvz_dy;
    FLOAT_P dvz_dz;
    FLOAT_P dd_vy_dy;
    FLOAT_P dd_vy_dz;
    FLOAT_P dd_vz_dy;
    FLOAT_P dd_vz_dz;
    FLOAT_P dd_vy_dydz;
    FLOAT_P dd_vz_dydz;
};

void calculate_common_central_derivatives(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct CommonCentralDerivatives *central_derivatives, int i, int j);

#endif // SPACIAL_DERIVATIVES_H__