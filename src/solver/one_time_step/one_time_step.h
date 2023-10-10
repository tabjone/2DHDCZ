#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "../functions.h"
#include "../rhs_functions/rhs_functions.h"
#include "shared_files.h"

void calculate_damping(FLOAT_P *damping_factor, struct GridInfo *grid_info, struct MpiInfo *mpi_info);
FLOAT_P get_dt(struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk1(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

#if DIMENSIONS == 1
FLOAT_P rk1_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
#elif DIMENSIONS == 2
FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
#elif DIMENSIONS == 3
FLOAT_P rk1_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk2_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
FLOAT_P rk3_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);
#endif

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep);


// OLD FUNCTIONS
FLOAT_P rk2(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk3(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);



#endif // ONE_TIME_STEP_H__