#include "one_time_step_3D.h"
#include <float.h>

FLOAT_P get_dt_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calulates the timestep by the CLF condition and the stability condition.

    Parameters
    ----------
    fg : struct
        A pointer to the ForegroundVariables2D struct.
    grid_info : struct
        A pointer to the GridInfo2D struct.
    dt_last : FLOAT_P
        The timestep used at the previous timestep.
    first_timestep : bool
        True if this is the first timestep, false otherwise.
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;

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
    #endif // TIME_ORDER, UPWIND_ORDER

    // CFL condition
    FLOAT_P dt_gridpoint; // dt at each grid point
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                dt_gridpoint = CFL_CUT * C_max / (fabs(fg->vx[i][j][k])/dx + fabs(fg->vy[i][j][k])/dy + fabs(fg->vz[i][j][k])/dz);
                if (dt_gridpoint < dt) 
                {
                    dt = dt_gridpoint;  // update dt if the new value is smaller
                }
            }
        }
    }
    
    // If this is the first timestep, we set dt to 0.1
    if (first_timestep && dt > 0.1 && MAX_DT > 0.1)
    {
        return 0.1;
    }

    // Stability condition, dt should be on the interval [0, 3/2]*dt_last
    if (dt > 1.5*dt_last && dt_last!=0.0)
    {
        dt = 1.5*dt_last;
    }
    
    // Maximum allowed timestep
    if (dt > MAX_DT)
    {
        dt = MAX_DT;
    }

    return dt;
}