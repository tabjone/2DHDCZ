#include <float.h>

void gauss_seidel_fast_second_order(double **b, double dz, double dx, int max_iterations, double tolerance)
{
    double a = 1.0/(dz*dz);
    double g = 1.0/(dx*dx);
    double c = 2.0*(a+g);

    // Initializing p and pnew to 0
    // When program runs properly: Set initial p to p1
    double pnew[nz][nx];
    double p[nz][nx];

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            pnew[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }

    int iter = 0;
    double max_difference, max_pnew;

    double abs_difference, abs_pnew;

    double tolerance_criteria = DBL_MAX;

    while (tolerance_criteria < tolerance)
    {
        max_difference = DBL_MAX;
        max_pnew = DBL_MAX;

        if (iter > max_iterations)
            {
                break;
            }

        // Copy p into pnew
        for (int i = 1; i < nz-1; i++)
        {
            for (j = 1; j < nx-1; j++)
            {
                p[i][j] = pnew[i][j]
            }
        }
        
        // Update pnew by eq. (54) from webpage (https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method) modified to dx!=dz
        for (int i = 1; i < nz-1; i++)
        {
            for (int j = 1; j < nx-1; j++)
            {
                pnew[i][j] = (a*(pnew[i+1][j]+p[i-1][j])+g*(p[i][j+1]+pnew[i][j-1])-b[i+nz_ghost][j])/c;

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
    #endif
}