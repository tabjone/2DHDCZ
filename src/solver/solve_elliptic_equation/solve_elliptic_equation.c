#include "solve_elliptic_equation.h"
#include "global_parameters.h"

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info)
{
    /*
    Solves the elliptic equation for the pressure field

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    fg_prev : struct ForegroundVariables
        Foreground variables at previous time step
    fg : struct ForegroundVariables
        Foreground variables at current time step
    grid_info : struct GridInfo
        Grid parameters
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;

    FLOAT_P **rhs;
    allocate_2D_array(&rhs, nz, ny);

    // Solve inside the grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            rhs[i][j] = rhs_elliptic_eq_2D(bg, fg_prev, grid_info, i+nz_ghost, j);
        }
    }

    // Solve boudaries
    for (int j = 0; j < ny; j++)
    { 
        // Bottom boundary
        rhs[0][j] = 0.0;
        // Top bundary
        rhs[nz-1][j] = 0.0;
    }

    gauss_seidel_2D(rhs, fg->p1, fg_prev->p1, grid_info);
    deallocate_2D_array(rhs);
}