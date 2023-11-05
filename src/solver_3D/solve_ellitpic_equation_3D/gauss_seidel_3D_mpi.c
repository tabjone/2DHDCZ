#include "solve_elliptic_equation_3D.h"
#include <float.h>

void gauss_seidel_3D_mpi(FLOAT_P ***b, FLOAT_P ***p1, FLOAT_P ***initial_p1, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
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
    allocate_3D_array(&pnew, nz+2, ny, nx);
    allocate_3D_array(&p, nz+2, ny, nx);
    
    // Initialize ghost cells
    for (int j = 0; j < ny; j++)
    {
        for (int k = 0; k < nx; k++)
        {
            pnew[0][j][k] = 0.0;
            pnew[nz+1][j][k] = 0.0;
            p[0][j][k] = 0.0;
            p[nz+1][j][k] = 0.0;
        }
    }
    
    // Initialize inside grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                pnew[i+1][j][k] = initial_p1[i+nz_ghost][j][k];
                p[i+1][j][k] = 0.0;
            }
        }
    }

    communicate_p_gauss_seidel_3D(pnew, grid_info, mpi_info);

    #if VERTICAL_BOUNDARY_TYPE == 2
        int i_start = 0;
        int i_end = nz;
    #else
        int i_start = 1;
        int i_end = nz+1;

        if (!mpi_info->has_neighbor_below) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    pnew[1][j][k] = 0.0;
                }
            }
        }
        if (!mpi_info->has_neighbor_above) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    pnew[nz][j][k] = 0.0;
                }
            }
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    // Tolerance parameters
    FLOAT_P max_difference, max_pnew;
    FLOAT_P abs_difference, abs_pnew;
    FLOAT_P tolerance_criteria = DBL_MAX;

    // Periodic boundary conditions
    int i_plus, i_minus, j_plus, j_minus, k_plus, k_minus;

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
        copy_3D_array(pnew, p, 0, nz+2, 0, ny, 0, nx);

        communicate_p_gauss_seidel_3D(p, grid_info, mpi_info);
    
        for (int i = i_start; i < i_end; i++)
        {
            #if VERTICAL_BOUNDARY_TYPE == 2
                i_plus = periodic_boundary(i+1, nz);
                i_minus = periodic_boundary(i-1, nz);
            #else
                i_plus = i+1;
                i_minus = i-1;
            #endif // VERTICAL_BOUNDARY_TYPE
            for (int j = 0; j < ny; j++)
            {
                j_plus = periodic_boundary(j+1, ny);
                j_minus = periodic_boundary(j-1, ny);
                for (int k = 0; k < nx; k++)
                {
                    k_plus = periodic_boundary(k+1, nx);
                    k_minus = periodic_boundary(k-1, nx);
                    // Gauss-Seidel
                    //pnew[i][j][k] = (b[i][j][k] - a*(p[i+1][j][k] + pnew[i-1][j][k]) - c*(p[i][j_plus][k] + pnew[i][j_minus][k]) - d*(p[i][j][k_plus] + pnew[i][j][k_minus]))/g;
                    // Jacobi
                    pnew[i][j][k] = (b[i][j][k] - a*(p[i_plus][j][k] + p[i_minus][j][k]) - c*(p[i][j_plus][k] + p[i][j_minus][k]) - d*(p[i][j][k_plus] + p[i][j][k_minus]))/g;
                
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
            for (int k = 0; k < nx; k++)
            {
                p1[i + nz_ghost][j][k] = pnew[i+1][j][k];
            }
        }
    }
    
    deallocate_3D_array(p);
    deallocate_3D_array(pnew);
}