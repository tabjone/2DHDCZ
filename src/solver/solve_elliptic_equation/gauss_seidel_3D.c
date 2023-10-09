#include "solve_elliptic_equation.h"
#include <float.h>

#if DIMENSIONS == 3
void gauss_seidel_3D(FLOAT_P ***b, FLOAT_P ***p1, FLOAT_P ***initial_p1, struct GridInfo *grid_info)
{
    /*
    Solves the elliptic equation for the pressure field using Gauss-Seidel method

    Parameters
    ----------
    b : FLOAT_P ***
        Right hand side of the elliptic equation
    p1 : FLOAT_P ***
        Pressure field at current time step
    initial_p1 : FLOAT_P ***
        Pressure field at previous time step
    grid_info : struct GridInfo
        Grid parameters
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Pre-calculating stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P c = 1.0/(dy*dy);
    FLOAT_P d = 1.0/(dx*dx);
    FLOAT_P g = -2.0*(a+c+d);

    // Initializing p and pnew
    FLOAT_P ***pnew, ***p;
    allocate_3D_array(&pnew, nz, ny, nx);
    allocate_3D_array(&p, nz, ny, nx);

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
            pnew[i][j][k] = initial_p1[i][j][k];
            p[i][j][k] = 0.0;
            }
        }
    }

    // Setting boundary conditions
    for (int j = 0; j < ny; j++)
    {
        for (int k = 0; k < nx; k++)
        {
            pnew[0][j][k] = 0.0;
            pnew[nz-1][j][k] = 0.0;
        }
    }

    // Tolerance parameters
    FLOAT_P max_difference, max_pnew;
    FLOAT_P abs_difference, abs_pnew;
    FLOAT_P tolerance_criteria = DBL_MAX;

    // Periodic boundary conditions
    int j_plus, j_minus, k_plus, k_minus;

    int iter = 0;
    while (tolerance_criteria > GS_TOL)
    {
        if (iter > GS_MAX_ITER)
        {
            break;
        }
        max_difference = 0.0;
        max_pnew = 0.0;

        // Copy pnew to p
        copy_3D_array(pnew, p, 0, nz, 0, ny, 0, nx);
        
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    j_plus = periodic_boundary(j+1, ny);
                    j_minus = periodic_boundary(j-1, ny);
                    k_plus = periodic_boundary(k+1, nx);
                    k_minus = periodic_boundary(k-1, nx);

                    pnew[i][j][k] = (b[i][j][k] - a*(p[i+1][j][k] + pnew[i-1][j][k]) - c*(p[i][j_plus][k] + pnew[i][j_minus][k]) - d*(p[i][j][k_plus] + pnew[i][j][k_minus]))/g;
                    
                    // Finding maximum absolute value of pnew
                    abs_pnew = fabs(pnew[i][j][k]);

                    if (abs_pnew > max_pnew)
                    {
                        max_pnew = abs_pnew;
                    }

                    // Finding maximum difference between p and pnew
                    abs_difference = fabs(pnew[i][j][k] - p[i][j][k]);
                    
                    if (abs_difference > max_difference)
                    {
                        max_difference = abs_difference;
                    }
                }
            }
        }
        tolerance_criteria = max_difference/max_pnew;
        iter++;
    }

    #if DEBUG == 1
        printf("Gauss-Seidel iterations: %d\n", iter);
    #endif // DEBUG

    // Updating p1
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                p1[i+nz_ghost][j][k] = pnew[i][j][k];
            }
        }
    }

    deallocate_3D_array(pnew);
    deallocate_3D_array(p);
}
#endif // DIMENSIONS == 3