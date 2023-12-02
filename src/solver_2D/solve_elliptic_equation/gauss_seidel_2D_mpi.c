#include "solve_elliptic_equation.h"
#include <float.h>

void gauss_seidel_2D_mpi(FLOAT_P **b, FLOAT_P **p1, FLOAT_P **initial_p1, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
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
    grid_info : struct GridInfo2D
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
    allocate_2D_array(&pnew, nz+2, ny);
    allocate_2D_array(&p, nz+2, ny);
    
    // Initialize ghost cells
    for (int j = 0; j < ny; j++)
    {
        pnew[0][j] = 0.0;
        pnew[nz+1][j] = 0.0;
        p[0][j] = 0.0;
        p[nz+1][j] = 0.0;
    }
    
    // Initialize inside grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            pnew[i+1][j] = initial_p1[i+nz_ghost][j];
            p[i+1][j] = 0.0;
        }
    }

    communicate_p_gauss_seidel(pnew, grid_info, mpi_info);

    int i_start = 1;
    int i_end = nz+1;

    // Non-periodic boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (!mpi_info->has_neighbor_below) // Boundary zero
        {
            i_start = 2;
            for (int j = 0; j < ny; j++)
            {
                pnew[1][j] = LOWER_PRESSURE_BOUNDARY;
            }
        }
        if (!mpi_info->has_neighbor_above) // Boundary zero
        {
            i_end = nz;
            for (int j = 0; j < ny; j++)
            {
                pnew[nz][j] = UPPER_PRESSURE_BOUNDARY;
            }
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
        copy_2D_array(pnew, p, 0, nz+2, 0, ny);

        communicate_p_gauss_seidel(p, grid_info, mpi_info);
    
        for (int i = i_start; i < i_end; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                j_plus = periodic_boundary(j+1, ny);
                j_minus = periodic_boundary(j-1, ny);
                // Gauss-Seidel
                //pnew[i][j] = (b[i][j] - a*(p[i_plus][j] + pnew[i_minus][j]) - c*(p[i][j_plus] + pnew[i][j_minus]))/g;
                // Jacobi
                //printf("i_minus=%d, i_plus=%d\n", i_minus, i_plus);
                pnew[i][j] = (b[i-1][j] - a*(p[i+1][j] + p[i-1][j]) - c*(p[i][j_plus] + p[i][j_minus]))/g;
                
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
        FLOAT_P reduced_max_difference;
        FLOAT_P reduced_max_pnew;

        MPI_Allreduce(&max_difference, &reduced_max_difference, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);
        // Find biggest max_pnew of all processes
        MPI_Allreduce(&max_pnew, &reduced_max_pnew, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);
        max_difference = reduced_max_difference;
        max_pnew = reduced_max_pnew;
        tolerance_criteria = max_difference/max_pnew;
        iter++;
    }

    #if DEBUG == 1
        if (mpi_info->rank == 0)
            printf("Gauss-Seidel iterations: %d\n", iter);
    #endif // DEBUG
    
    // Updating p1
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            p1[i + nz_ghost][j] = pnew[i+1][j];
        }
    }
    
    deallocate_2D_array(p);
    deallocate_2D_array(pnew);
}