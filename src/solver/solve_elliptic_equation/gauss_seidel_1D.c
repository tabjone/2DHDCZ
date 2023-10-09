#include "solve_elliptic_equation.h"
#include <float.h>

#if DIMENSIONS == 1
void gauss_seidel_1D(FLOAT_P *b, FLOAT_P *p1, FLOAT_P *initial_p1, struct GridInfo *grid_info)
{
    /*
    Solves the elliptic equation for the pressure field using Gauss-Seidel method

    Parameters
    ----------
    b : FLOAT_P *
        Right hand side of the elliptic equation
    p1 : FLOAT_P *
        Pressure field at current time step
    initial_p1 : FLOAT_P *
        Pressure field at previous time step
    grid_info : struct GridInfo
        Grid parameters
    */
    
    // Getting grid info
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dz = grid_info->dz;

    // Pre-calculating stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P c = -2.0*a;

    // Initializing p and pnew
    FLOAT_P pnew[nz];
    FLOAT_P p[nz];

    for (int i = 0; i < nz; i++)
    {
        pnew[i] = initial_p1[i];
        p[i] = 0.0;
    }

    // Setting boundary conditions
    pnew[0] = 0.0;
    pnew[nz-1] = 0.0;

    // Tolerance parameters
    FLOAT_P max_difference, max_pnew;
    FLOAT_P abs_difference, abs_pnew;
    FLOAT_P tolerance_criteria = DBL_MAX;

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
        copy_1D_array(pnew, p, 0, nz);

        for (int i = 1; i < nz-1; i++)
        {
            pnew[i] = (b[i] - a*(p[i+1]+pnew[i-1]))/c;

            // Finding maximum absolute value of pnew
            abs_pnew = fabs(pnew[i]);
            
            if (abs_pnew > max_pnew)
            {
                max_pnew = abs_pnew;
            }

            // Finding maximum difference between p and pnew
            abs_difference = fabs(pnew[i] - p[i]);

            if (abs_difference > max_difference)
            {
                max_difference = abs_difference;
            }
        }
        tolerance_criteria = max_difference/max_pnew;
        iter++;
    }

    // Updating p1
    for (int i = 0; i < nz; i++)
    {
        p1[i+nz_ghost] = pnew[i];
    }

    #if DEBUG == 1
        printf("Gauss-Seidel iterations: %d\n", iter);
    #endif // DEBUG
}
#endif // DIMENSIONS == 1