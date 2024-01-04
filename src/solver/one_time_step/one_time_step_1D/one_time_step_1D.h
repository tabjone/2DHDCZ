#ifndef ONE_TIME_STEP_1D_H__
#define ONE_TIME_STEP_1D_H__

#include <stdbool.h>
#include "global_float_precision.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "iterative_solver_module/iterative_solver_1D/iterative_solver_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"
#include "MPI_module/MPI_module.h"
#include "solver/rhs_functions/rhs_functions_1D/rhs_functions_1D.h"

FLOAT_P get_dt_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk1_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc, FLOAT_P dt_last, bool first_timestep);

FLOAT_P one_time_step_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc, FLOAT_P dt_last, bool first_timestep);

void solve_elliptic_equation_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc);

#endif // ONE_TIME_STEP_1D_H__