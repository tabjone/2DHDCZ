#ifndef ONE_TIME_STEP_3D_H__
#define ONE_TIME_STEP_3D_H__

#include "global_parameters.h"
#include "shared_files.h"
#include "../equations_3D/equations_3D.h"
#include "../rhs_functions_3D/rhs_functions_3D.h"
#include "../solve_ellitpic_equation_3D/solve_elliptic_equation_3D.h"
#include "../boundary_3D/boundary_3D.h"

FLOAT_P get_dt_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P one_time_step_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk3_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);

#endif // ONE_TIME_STEP_3D_H__