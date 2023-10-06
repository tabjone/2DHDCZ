#include <float.h>
#include "solve_elliptic_equation.h"

void gauss_seidel_fast_second_order(FLOAT_P **b, FLOAT_P **p1, FLOAT_P **initial_p1, struct GridInfo *grid_info)
{
    int nx = grid_info->nx;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P g = 1.0/(dx*dx);
    FLOAT_P c = 2.0*(a+g);

    // Initializing p and pnew to 0
    // When program runs properly: Set initial p to p1
    FLOAT_P pnew[nz][nx];
    FLOAT_P p[nz][nx];

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            pnew[i][j] = initial_p1[i+nz_ghost][j];
            p[i][j] = 0.0;
        }
    }

    for (int j = 0; j < nx; j++)
    {
        pnew[0][j] = 0.0;
        pnew[nz-1][j] = 0.0;
    }

    int iter = 0;
    FLOAT_P max_difference, max_pnew;

    FLOAT_P abs_difference, abs_pnew;

    FLOAT_P tolerance_criteria = DBL_MAX;

    int j_plus, j_minus;

    while (tolerance_criteria > GS_TOL)
    {
        max_difference = 0.0;
        max_pnew = 0.0;

        if (iter > GS_MAX_ITER)
            {
                break;
            }

        // Copy pnew to p
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                p[i][j] = pnew[i][j];
            }
        }
        
        // Update pnew by eq. (54) from webpage (https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method) modified to dx!=dz
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                j_plus = periodic_boundary(j+1, nx);
                j_minus = periodic_boundary(j-1, nx);

                pnew[i][j] = (a*(pnew[i+1][j]+p[i-1][j])+g*(p[i][j_plus]+pnew[i][j_minus])-b[i+nz_ghost][j])/c;
                //pnew[i][j] = (a*(p[i+1][j]+p[i-1][j])+g*(p[i][j_plus]+p[i][j_minus])-b[i+nz_ghost][j])/c;


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

    // Updating p1
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            p1[i+nz_ghost][j] = pnew[i][j];
        }
    }

    #if DEBUG == 1
    printf("Gauss-Seidel iterations: %d\n", iter);
    #endif
}