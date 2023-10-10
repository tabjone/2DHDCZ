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

    #if DIMENSIONS == 1
        // Getting grid info
        int nz = grid_info->nz;
        int nz_ghost = grid_info->nz_ghost;

        FLOAT_P *rhs;
        allocate_1D_array(&rhs, nz);

        // Solve inside the grid
        for (int i = 1; i < nz-1; i++)
        {
            rhs[i] = rhs_elliptic_eq_1D(bg, fg_prev, grid_info, i+nz_ghost);
        }

        // Solve boudaries
        // Bottom boundary
        rhs[0] = 0.0;
        // Top bundary
        rhs[nz-1] = 0.0;

        gauss_seidel_1D(rhs, fg->p1, fg_prev->p1, grid_info);
        deallocate_1D_array(rhs);

    #elif DIMENSIONS == 2
        // Getting grid info
        int ny = grid_info->ny;
        int nz_ghost = grid_info->nz_ghost;
        int nz = grid_info->nz;

        FLOAT_P **rhs;
        allocate_2D_array(&rhs, nz, ny);

        // Solve inside the grid
        for (int i = 1; i < nz-1; i++)
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

    #elif DIMENSIONS == 3
        // Getting grid info
        int nx = grid_info->nx;
        int ny = grid_info->ny;
        int nz_ghost = grid_info->nz_ghost;
        int nz = grid_info->nz;

        FLOAT_P ***rhs;
        allocate_3D_array(&rhs, nz, ny, nx);

        // Solve inside the grid
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    rhs[i][j][k] = rhs_elliptic_eq_3D(bg, fg_prev, grid_info, i+nz_ghost, j, k);
                }
            }
        }

        // Solve boudaries
        for (int j = 0; j < ny; j++)
        { 
            for (int k = 0; k < nx; k++)
            {
                // Bottom boundary
                rhs[0][j][k] = 0.0;
                // Top bundary
                rhs[nz-1][j][k] = 0.0;
            }
        }

        gauss_seidel_3D(rhs, fg->p1, fg_prev->p1, grid_info);
        deallocate_3D_array(rhs);
    #endif // DIMENSIONS
}