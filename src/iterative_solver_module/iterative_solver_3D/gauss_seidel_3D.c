#include <mpi.h>
#include <float.h>
#include "global_float_precision.h"
#include "global_boundary.h"
#include "global_parameters.h"
#include "MPI_module/MPI_module.h"
#include "array_utilities/array_copy/array_copy.h"
#include <math.h>

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

void gauss_seidel_3D(FLOAT_P ***rhs, FLOAT_P ***current_solution, FLOAT_P ***previous_solution, int nz, int nz_ghost, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, struct MpiInfo *mpi_info)
{
    /*
    Solves the Poisson equation using the Gauss-Seidel method

    Parameters
    ----------
    rhs : FLOAT_P **
        Right hand side of the Poisson equation
    current_solution : FLOAT_P **
        Solution of the Poisson equation at current iteration
    previous_solution : FLOAT_P **
        Solution of the Poisson equation at previous iteration
    nz : int
        Number of grid points in z-direction
    nz_ghost : int
        Number of ghost points in z-direction
    ny : int
        Number of grid points in y-direction
    nx : int
        Number of grid points in x-direction
    dz : FLOAT_P
        Grid spacing in z-direction
    dy : FLOAT_P
        Grid spacing in y-direction
    dx : FLOAT_P
        Grid spacing in x-direction
    mpi_info : struct MpiInfo
        Struct containing MPI information
    */

    // Stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P c = 1.0/(dy*dy);
    FLOAT_P d = 1.0/(dx*dx);
    FLOAT_P g = -2.0*(a+c+d);

    // Tolerance parameters
    FLOAT_P local_abs_difference, my_global_abs_difference, global_abs_difference;;
    FLOAT_P local_abs_current, my_global_abs_current, global_abs_current;
    FLOAT_P tolerance_criteria = DBL_MAX;

    // Periodic boundary conditions in y-direction
    int j_plus, j_minus;
    int k_plus, k_minus;

    int i_start = 1; // For periodic boundary we also solve on the boundary
    int i_end = nz+1; // For periodic boundary we also solve on the boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (!mpi_info->has_neighbor_below)
        {i_start = 2;} // When not periodic, we do not solve on the boundary
        if (!mpi_info->has_neighbor_above)
        {i_end = nz;} // When not periodic, we do not solve on the boundary
    #endif // VERTICAL_BOUNDARY_TYPE

    int iter = 0;
    while (tolerance_criteria > ITERATIVE_SOLVER_TOLERANCE)
    {
        if (iter > ITERATIVE_SOLVER_MAX_ITERATIONS)
        {
            break;
        }
        my_global_abs_difference = 0.0;
        my_global_abs_current = 0.0;

        // Copy current_solution to previous_solution and communicate ghost cells
        copy_3D_array(current_solution, previous_solution, 0, nz+2, 0, ny, 0, nx);
        communicate_3D_ghost_above_below(previous_solution, mpi_info, nz, 1, ny, nx);
    
        for (int i = i_start; i < i_end; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                j_plus = periodic_boundary(j+1, ny);
                j_minus = periodic_boundary(j-1, ny);
                for (int k = 0; k < nx; k++)
                {
                    k_plus = periodic_boundary(k+1, nx);
                    k_minus = periodic_boundary(k-1, nx);

                    current_solution[i][j][k] = 
                    (rhs[i-1][j][k] 
                    -a*(previous_solution[i+1][j][k] + current_solution[i-1][j][k]) 
                    -c*(previous_solution[i][j_plus][k] + current_solution[i][j_minus][k])
                    -d*(previous_solution[i][j][k_plus] + current_solution[i][j][k_minus]))/g;
                
                    // Finding maximum absolute value of current_solution
                    local_abs_current = fabs(current_solution[i][j][k]);

                    if (local_abs_current > my_global_abs_current)
                    {
                        my_global_abs_current = local_abs_current;
                    }

                    // Finding maximum difference between previous_solution and current_solution
                    local_abs_difference = fabs(current_solution[i][j][k] - previous_solution[i][j][k]);
                    
                    if (local_abs_difference > my_global_abs_difference)
                    {
                        my_global_abs_difference = local_abs_difference;
                    }
                }
            }
        }

        // Find biggest abs_difference and abs_current of all processes
        MPI_Allreduce(&my_global_abs_difference, &global_abs_difference, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&my_global_abs_current, &global_abs_current, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);

        tolerance_criteria = global_abs_difference/global_abs_current;
        iter++;
    }

    if (mpi_info->rank == 0)
        {printf("Gauss-Seidel converged after iterations: %d\n", iter);}
}