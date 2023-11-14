#include "solve_elliptic_equation.h"
#include "global_parameters.h"

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc)
{
    /*
    Solves the elliptic equation for the pressure field

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    fg_prev : struct ForegroundVariables2D
        Foreground variables at previous time step
    fg : struct ForegroundVariables2D
        Foreground variables at current time step
    grid_info : struct GridInfo2D
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
            rhs[i][j] = rhs_elliptic_eq_2D(bg, fg_prev, grid_info, precalc, i+nz_ghost, j);
        }
    }
    #if MPI_ON == 0
        gauss_seidel_2D(rhs, fg->p1, fg_prev->p1, grid_info);
    #elif MPI_ON == 1
        gauss_seidel_2D_mpi(rhs, fg->p1, fg_prev->p1, grid_info, mpi_info);
    #endif // MPI_ON
    deallocate_2D_array(rhs);
}