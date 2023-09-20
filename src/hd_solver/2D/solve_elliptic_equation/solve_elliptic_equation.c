#include "solve_elliptic_equation.h"

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg)
{
    int nx = fg->nx;
    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;
    int nz = fg->nz;

    double dx = fg->dx;
    double dz = fg->dz;

    double **rhs;
    // CHECK THIS: THIS CAN BE OF SIZE nz, nx INSTEAD OF nz_full, nx !!!
    allocate_2D_array(&rhs, nz_full, nx);

    // Solve inside the grid
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            rhs[i][j] = - rhs_elliptic_eq(bg, fg, i, j);
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

    gauss_seidel_fast_second_order(rhs, dz, dx, nx, nz, nz_ghost, 1e5, 1e-6, fg->p1);
    deallocate_2D_array(rhs);

    /*
    // Create rhs array without ghost cells
    double *rhs_without_ghost_cells;
    allocate_1D_array(&rhs_without_ghost_cells, nx*nz);
    for (int i = 0; i < nx*nz; i++)
    {
        rhs_without_ghost_cells[i] = rhs[i/nx+nz_ghost][i%nx];
    }
    // Create A matrix
    double **A;

    int N = nx*nz;

    allocate_2D_array(&A, N, N);
    int diagonal_pos = 0;
    //printf("N = %d\n", N);
    //printf("nx*nz = %d\n", nx*nz);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = 0.0;
        }
        A[i][diagonal_pos] = -2.0*(1.0/(dx*dx)+1.0/(dz*dz));
        if (i>0)
        { A[i][diagonal_pos-1] = 1.0/(dx*dx); }
        if (i<nx*nz-1)
        { A[i][diagonal_pos+1] = 1.0/(dx*dx); }
        if (i>nx-1)
        { A[i][diagonal_pos-nx] = 1.0/(dz*dz); }
        if (i<nx*nz-nx)
        { A[i][diagonal_pos+nx] = 1.0/(dz*dz); }
        diagonal_pos += 1;
    }

    double *x;
    allocate_1D_array(&x, N);
    gauss_seidel(A, rhs_without_ghost_cells, x, N, 100, 1e-6);
 
    // Put x into p1
    for (int i = 0; i < N; i++)
    {
        fg->p1[i/nx+nz_ghost][i%nx] = x[i];
    }

    deallocate_1D_array(rhs_without_ghost_cells);
    deallocate_2D_array(rhs);
    deallocate_2D_array(A);
    deallocate_1D_array(x);
    */
}