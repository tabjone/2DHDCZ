#include <mpi.h>
#include <float.h>
#include "global_float_precision.h"
#include "global_boundary.h"
#include "global_parameters.h"
#include "MPI_module/MPI_module.h"
#include "array_utilities/array_copy/array_copy.h"
#include <math.h>

void jacobi_1D(FLOAT_P *rhs, FLOAT_P *current_solution, FLOAT_P *previous_solution, int nz, int nz_ghost, FLOAT_P dz, struct MpiInfo *mpi_info)
{
    /*
    Solves the Poisson equation using the Jacobi method

    Parameters
    ----------
    rhs : FLOAT_P *
        Right hand side of the Poisson equation
    current_solution : FLOAT_P *
        Solution of the Poisson equation at current iteration
    previous_solution : FLOAT_P *
        Solution of the Poisson equation at previous iteration
    nz : int
        Number of grid points in z-direction
    nz_ghost : int
        Number of ghost points in z-direction
    dz : FLOAT_P
        Grid spacing in z-direction
    mpi_info : struct MpiInfo
        Struct containing MPI information
    */

    // Stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P g = -2.0*(a);

    // Tolerance parameters
    FLOAT_P local_abs_difference, my_global_abs_difference, global_abs_difference;;
    FLOAT_P local_abs_current, my_global_abs_current, global_abs_current;
    FLOAT_P tolerance_criteria = DBL_MAX;

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
        copy_1D_array(current_solution, previous_solution, 0, nz+2);
        communicate_1D_ghost_above_below(previous_solution, mpi_info, nz, 1);
    
        for (int i = i_start; i < i_end; i++)
        {
            current_solution[i] = (rhs[i-1] - a*(previous_solution[i+1] + previous_solution[i-1]))/g;
            
            // Finding maximum absolute value of current_solution
            local_abs_current = fabs(current_solution[i]);

            if (local_abs_current > my_global_abs_current)
            {
                my_global_abs_current = local_abs_current;
            }

            // Finding maximum difference between previous_solution and current_solution
            local_abs_difference = fabs(current_solution[i] - previous_solution[i]);
            
            if (local_abs_difference > my_global_abs_difference)
            {
                my_global_abs_difference = local_abs_difference;
            }
        }

        // Find biggest abs_difference and abs_current of all processes
        MPI_Allreduce(&my_global_abs_difference, &global_abs_difference, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&my_global_abs_current, &global_abs_current, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);

        tolerance_criteria = global_abs_difference/global_abs_current;
        iter++;
    }

    if (mpi_info->rank == 0)
        {printf("Jacobi converged after iterations: %d\n", iter);}
}