#ifndef ONE_TIME_STEP_3D_H__
#define ONE_TIME_STEP_3D_H__

#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_3D/iterative_solver_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_struct_3D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_3D/rhs_functions_3D.h"

FLOAT_P get_dt_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk1_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc, FLOAT_P dt_last, bool first_timestep);

FLOAT_P one_time_step_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc, FLOAT_P dt_last, bool first_timestep);

void solve_elliptic_equation_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc);

#endif // ONE_TIME_STEP_3D_H__