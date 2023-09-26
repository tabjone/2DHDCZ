#include "solve_elliptic_equation.h"

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info)
{
    int nx = grid_info->nx;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    FLOAT_P **rhs;
    // CHECK THIS: THIS CAN BE OF SIZE nz, nx INSTEAD OF nz_full, nx !!!
    allocate_2D_array(&rhs, nz_full, nx);

    // Solve inside the grid
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            rhs[i][j] = - rhs_elliptic_eq(bg, fg_prev, grid_info, i, j);
        }
    }
    

    // Solve boudaries
    for (int j = 0; j < nx; j++)
    { 
        // Bottom boundary
        rhs[nz_ghost][j] = 0.0;
        // Top bundary
        rhs[nz_full-nz_ghost-1][j] = 0.0;
    }

    gauss_seidel_fast_second_order(rhs, fg->p1, fg_prev->p1, grid_info);
    deallocate_2D_array(rhs);
}