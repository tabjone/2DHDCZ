#include <mpi.h>
#include <float.h>
#include "global_float_precision.h"
#include "global_boundary.h"
#include "global_parameters.h"
#include "MPI_module/MPI_module.h"
#include "array_utilities/array_copy/array_copy.h"
#include <math.h>

#include "periodic_boundary.h"

void jacobi_2D_spherical(FLOAT_P **rhs, FLOAT_P **current_solution, FLOAT_P **previous_solution, FLOAT_P *r_array, int nz, int nz_ghost, int ny, FLOAT_P dz, FLOAT_P dy, struct MpiInfo *mpi_info)
{
    /*
    Solves the Poisson equation using the Jacobi method

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
    dz : FLOAT_P
        Grid spacing in z-direction
    dy : FLOAT_P
        Grid spacing in y-direction
    mpi_info : struct MpiInfo
        Struct containing MPI information
    */

    // Stencil values
    FLOAT_P a = 1.0/(dz*dz);
    FLOAT_P b = 1.0/(2.0*dz);
    FLOAT_P c = 1.0/(dy*dy);
    FLOAT_P d;

    FLOAT_P r; // Radius

    // Tolerance parameters
    FLOAT_P iterative_difference;
    FLOAT_P global_iterative_difference_sum_squares, local_iterative_difference_sum_squares;
    FLOAT_P global_current_solution_sum_squares, local_current_solution_sum_squares;
    FLOAT_P current_solution_norm, iterative_difference_norm;

    FLOAT_P tolerance_criteria = DBL_MAX;

    // Periodic boundary conditions in y-direction
    int j_plus, j_minus;

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
        local_iterative_difference_sum_squares = 0.0;
        local_current_solution_sum_squares = 0.0;

        // Copy current_solution to previous_solution and communicate ghost cells
        copy_2D_array(current_solution, previous_solution, 0, nz+2, 0, ny);
        communicate_2D_ghost_above_below(previous_solution, mpi_info, nz, 1, ny);
    
        for (int i = i_start; i < i_end; i++)
        {
            r = r_array[i];
            d = 2.0*a + 2.0*c/(r*r);
            for (int j = 0; j < ny; j++)
            {
                j_plus = periodic_boundary(j+1, ny);
                j_minus = periodic_boundary(j-1, ny);

                current_solution[i][j] = (-rhs[i-1][j] 
                + a*(previous_solution[i+1][j] + previous_solution[i-1][j]) 
                + b/r*(previous_solution[i+1][j] - previous_solution[i-1][j]) 
                + c/(r*r)*(previous_solution[i][j_plus] + previous_solution[i][j_minus]))/d;
                
                iterative_difference = current_solution[i][j] - previous_solution[i][j];
                local_iterative_difference_sum_squares += iterative_difference*iterative_difference;
                local_current_solution_sum_squares += current_solution[i][j]*current_solution[i][j];
            }
        }
        // Global sum of squares of difference between previous_solution and current_solution
        MPI_Allreduce(&local_iterative_difference_sum_squares, &global_iterative_difference_sum_squares, 1, MPI_FLOAT_P, MPI_SUM, MPI_COMM_WORLD);
        // Global sum of squares of current_solution
        MPI_Allreduce(&local_current_solution_sum_squares, &global_current_solution_sum_squares, 1, MPI_FLOAT_P, MPI_SUM, MPI_COMM_WORLD);

        iterative_difference_norm = sqrt(global_iterative_difference_sum_squares);
        current_solution_norm = sqrt(global_current_solution_sum_squares);

        tolerance_criteria = iterative_difference_norm / current_solution_norm;
        iter++;
    }

    if (mpi_info->rank == 0)
        {printf("Jacobi converged after iterations: %d\n", iter);}
}