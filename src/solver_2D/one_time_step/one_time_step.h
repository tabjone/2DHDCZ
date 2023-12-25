#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "global_float_precision.h"
#include "../../data_structures/background_data/background_variables_struct.h"
#include "../../data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "../../data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include <stdbool.h>

FLOAT_P get_dt(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep);

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc)

#endif // ONE_TIME_STEP_H__