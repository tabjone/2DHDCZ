#include "solve_elliptic_equation_3D.h"

void solve_elliptic_equation_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
{
    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;

    FLOAT_P ***rhs;
    allocate_3D_array(&rhs, nz, ny, nx);

    // Solve inside the grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                rhs[i][j][k] = rhs_elliptic_eq_3D(bg, fg_prev, grid_info, i+nz_ghost, j, k);
            }
        }
    }

    #if MPI_ON == 0
        gauss_seidel_3D(rhs, fg->p1, fg_prev->p1, grid_info);
    #elif MPI_ON == 1
        gauss_seidel_3D_mpi(rhs, fg->p1, fg_prev->p1, grid_info, mpi_info);
    #endif // MPI_ON

    deallocate_3D_array(rhs);

}