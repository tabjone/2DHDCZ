#include "solve_elliptic_equation.h"

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg)
{
    int nx = fg->nx;
    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;

    double **rhs; // Maybe I should not include ghost cells on this. Will make it easier later because this is b.
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
        // Top boundary
        rhs[nz_ghost][j] = rhs_elliptic_eq_vertical_boundary(bg, fg, nz_ghost, j);
        // Bottom boundary
        rhs[nz_full-nz_ghost-1][j] = rhs_elliptic_eq_vertical_boundary(bg, fg, nz_full-nz_ghost-1, j);

        for (int k = 0; k < nz_ghost; k++)
        {
            rhs[k][j] = rhs[nz_ghost][j];
            rhs[nz_full-k-1][j] = rhs[nz_full-nz_ghost-1][j];
        }
    }
    /*
    double *x, *b;
    double **A;
    int maxIterations = 1000;
    double tolerance = 1e-6;

    int nz = fg->nz;
    
    allocate_1D_array(&x, nx*nz);
    allocate_1D_array(&b, nx*nz);
    allocate_2D_array(&A, nx*nz, nx*nz);

    // Fill x
    for (int i = 0; i < nx*nz; i++)
    {
        x[i] = fg->p1[i/nx+nz_ghost][i%nx];
        b[i] = rhs[i/nx+nz_ghost][i%nx];
    }

    // Fill A
    int z_pos = 0;
    for (int i = 0; i < nx*nz; i++)
    {
        for (int j = 0; j < nx*nz; j++)
        {
            A[i][j] = 0.0;
        }
        A[i][z_pos] = -2.0*(1.0/(dx*dx)+1.0/(dz*dz));
        if (i>0)
        { A[i][z_pos-1] = 1.0/(dx*dx); }
        if (i<nx*nz-1)
        { A[i][z_pos+1] = 1.0/(dx*dx); }
        if (i>nx-1)
        { A[i][z_pos-nx] = 1.0/(dz*dz); }
        if (i<nx*nz-nx)
        { A[i][z_pos+nx] = 1.0/(dz*dz); }
        z_pos += 1;
    }

    // Print matrix
    for (int i = 0; i < nx*nz; i++)
    {
        for (int j = 0; j < nx*nz; j++)
        {
            //printf("%f ", A[i][j]);
        }
        //printf("\n");
    }


    
    // Now use Gauss-Seidel to solve the equation
    gaussSeidel(A, b, x, nx*nz, maxIterations, tolerance);
    
    deallocate_1D_array(b);
    deallocate_2D_array(A);
    deallocate_1D_array(x);*/
    deallocate_2D_array(rhs);
}