#ifndef ONE_TIME_STEP_2D_H__
#define ONE_TIME_STEP_2D_H__

#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_2D/iterative_solver_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_2D/rhs_functions_2D.h"

FLOAT_P get_dt_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);

FLOAT_P one_time_step_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);

void solve_elliptic_equation_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc);

#endif // ONE_TIME_STEP_2D_H__