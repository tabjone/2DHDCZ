#include "one_time_step.h"
#include <float.h>

FLOAT_P get_dt(struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dx = grid_info->dx;

    FLOAT_P dt;
    
    // We set dt to the maximum representable finite floating-point until we find a smaller value
    #if FLOAT_PRECISION == 0
        dt = FLT_MAX;
    #elif FLOAT_PRECISION == 1
        dt = DBL_MAX;
    #elif FLOAT_PRECISION == 2
        dt = LDBL_MAX;
    #endif

    FLOAT_P C_max; // CFL number fround from Von Neumann stability analysis

    #if TIME_ORDER == 1 && UPWIND_ORDER == 1
        C_max = 1.0;
    #elif TIME_ORDER == 1 && UPWIND_ORDER == 2
        #error "UNSTABLE"
    #elif TIME_ORDER == 2 && UPWIND_ORDER == 1
        C_max = 1.0;
    #elif TIME_ORDER == 2 && UPWIND_ORDER == 2
        C_max = 0.5;
    #elif TIME_ORDER == 3 && UPWIND_ORDER == 1
        C_max = 1.25;
    #elif TIME_ORDER == 3 && UPWIND_ORDER == 2
        C_max = 0.5;
    #else
        #error "TIME_ORDER must be either 1, 2, or 3, UPWIND_ORDER must be either 1 or 2"
    #endif

    // CFL condition
    FLOAT_P dt_gridpoint; // dt at each grid point
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dt_gridpoint = CFL_CUT * C_max / (fabs(fg->vx[i][j])/dx + fabs(fg->vz[i][j])/dz);
            if (dt_gridpoint < dt) 
            {
                dt = dt_gridpoint;  // update dt if the new value is smaller
            }
        }
    }
    // If this is the first timestep, we set dt to 0.1
    if (first_timestep && dt > 0.1)
    {
        return 0.1;
    }

    // Stability condition, dt should be on the interval [1/2, 2]*dt_last
    if (dt > 1.5*dt_last)
    {
        dt = 1.5*dt_last;
    }

    if (dt > MAX_DT)
    {
        dt = MAX_DT;
    }

    


    return dt;
}