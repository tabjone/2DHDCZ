#include "global_parameters.h"
#include "global_float_precision.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "solver/one_time_step/one_time_step_2D/one_time_step_2D.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

void step_entropy_velocity_2D(struct ForegroundVariables2D *fg, struct ForegroundVariables2D *fg_prev, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, FLOAT_P rhs_s1, FLOAT_P rhs_vy, FLOAT_P rhs_vz, FLOAT_P dt, int i, int j)
{
    int ny = grid_info->ny;
    int j_plus = periodic_boundary(j+1, ny);
    int j_minus = periodic_boundary(j-1, ny);

    fg->s1[i][j] = dt * rhs_s1;
    fg->vy[i][j] = dt * rhs_vy;
    fg->vz[i][j] = dt * rhs_vz;

    #if REMOVE_AVG_VZ_X == 1
        fg->vz[i][j] -= precalc->vz_horizontal_average[i]/AVG_VZ_FACTOR * dt;
    #endif

    #if LAX_S1 == 1
        fg->s1[i][j] += 1.0/(LAX_PARAM_S1+4.0)*(LAX_PARAM_S1*fg_prev->s1[i][j] + fg_prev->s1[i+1][j] + fg_prev->s1[i-1][j] + fg_prev->s1[i][j_plus] + fg_prev->s1[i][j_minus]);
    #else
        fg->s1[i][j] += fg_prev->s1[i][j];
    #endif

    #if LAX_VY == 1
        fg->vy[i][j] += 1.0/(LAX_PARAM_VY+4.0)*(LAX_PARAM_VY*fg_prev->vy[i][j] + fg_prev->vy[i+1][j] + fg_prev->vy[i-1][j] + fg_prev->vy[i][j_plus] + fg_prev->vy[i][j_minus]);
    #else
        fg->vy[i][j] += fg_prev->vy[i][j];
    #endif

    #if LAX_VZ == 1
        fg->vz[i][j] += 1.0/(LAX_PARAM_VZ+4.0)*(LAX_PARAM_VZ*fg_prev->vz[i][j] + fg_prev->vz[i+1][j] + fg_prev->vz[i-1][j] + fg_prev->vz[i][j_plus] + fg_prev->vz[i][j_minus]);
    #else
        fg->vz[i][j] += fg_prev->vz[i][j];
    #endif
}