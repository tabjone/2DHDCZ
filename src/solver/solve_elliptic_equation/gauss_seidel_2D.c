#include "solve_elliptic_equation.h"
#include <float.h>

#if DIMENSIONS == 2
void gauss_seidel_2D(FLOAT_P **b, FLOAT_P **p1, FLOAT_P **initial_p1, struct GridInfo *grid_info)
{
    /*
    Solves the elliptic equation for the pressure field using Gauss-Seidel method

    Parameters
    ----------
    b : FLOAT_P **
        Right hand side of the elliptic equation
    p1 : FLOAT_P **
        Pressure field at current time step
    initial_p1 : FLOAT_P **
        Pressure field at previous time step
    grid_info : struct GridInfo
        Grid parameters
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Pre-calculating stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P c = 1.0/(dy*dy);
    FLOAT_P g = -2.0*(a+c);

    // Initializing p and pnew
    FLOAT_P **pnew, **p;
    allocate_2D_array(&pnew, nz, ny);
    allocate_2D_array(&p, nz, ny);

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            pnew[i][j] = initial_p1[i][j];
            p[i][j] = 0.0;
        }
    }

    #if VERTICAL_BOUNDARY_TYPE != 2
        // Setting boundary conditions
        for (int j = 0; j < ny; j++)
        {
            pnew[0][j] = 0.0;
            pnew[nz-1][j] = 0.0;
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    // Tolerance parameters
    FLOAT_P max_difference, max_pnew;
    FLOAT_P abs_difference, abs_pnew;
    FLOAT_P tolerance_criteria = DBL_MAX;

    // Periodic boundary conditions
    int j_plus, j_minus;

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
        copy_2D_array(pnew, p, 0, nz, 0, ny);
        
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                j_plus = periodic_boundary(j+1, ny);
                j_minus = periodic_boundary(j-1, ny);

                pnew[i][j] = (b[i][j] - a*(p[i+1][j] + pnew[i-1][j]) - c*(p[i][j_plus] + pnew[i][j_minus]))/g;
                
                // Finding maximum absolute value of pnew
                abs_pnew = fabs(pnew[i][j]);

                if (abs_pnew > max_pnew)
                {
                    max_pnew = abs_pnew;
                }

                // Finding maximum difference between p and pnew
                abs_difference = fabs(pnew[i][j] - p[i][j]);
                
                if (abs_difference > max_difference)
                {
                    max_difference = abs_difference;
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
            p1[i + nz_ghost][j] = pnew[i][j];
        }
    }

    deallocate_2D_array(p);
    deallocate_2D_array(pnew);
}
#endif // DIMENSIONS == 2