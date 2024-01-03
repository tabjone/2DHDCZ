#include "global_float_precision.h"
#include "global_boundary.h"
#include "global_parameters.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "MPI_module/MPI_module.h"
#include "iterative_solver_1D.h"

void iterative_solver_1D(FLOAT_P *rhs, FLOAT_P *final_solution, FLOAT_P *initial_guess, int nz, int nz_ghost, FLOAT_P dz, struct MpiInfo *mpi_info)
{
    /*
    Solves the Poisson equation using an iterative method (Jacobi, Gauss-Seidel)

    Parameters
    ----------
    rhs : FLOAT_P **
        Right hand side of the Poisson equation
    final_solution : FLOAT_P **
        Array where solution will be stored
    initial_guess : FLOAT_P **
        Initial guess for the solution
    nz : int
        Number of grid points in z-direction
    nz_ghost : int
        Number of ghost points in z-direction
    dz : FLOAT_P
        Grid spacing in z-direction
    mpi_info : struct MpiInfo
        Struct containing MPI information
    */

    // Initializing previous_solution and current_solution
    FLOAT_P **current_solution, **previous_solution;
    allocate_1D_array(&current_solution, nz+2);
    allocate_1D_array(&previous_solution, nz+2);
    
    // Initialize ghost cells
    current_solution[0] = 0.0;
    current_solution[nz+1] = 0.0;
    
    // Initialize inside grid
    for (int i = 0; i < nz; i++)
    {
        current_solution[i+1] = initial_guess[i+nz_ghost];
    }

    communicate_1D_ghost_above_below(current_solution, mpi_info, nz, 1);

    // Non-periodic boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (!mpi_info->has_neighbor_below) // Boundary zero
        {
            current_solution[1] = LOWER_PRESSURE_BOUNDARY;
        }
        if (!mpi_info->has_neighbor_above) // Boundary zero
        {
            current_solution[nz] = UPPER_PRESSURE_BOUNDARY;
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if ITERATIVE_SOLVER_TYPE == 0
        jacobi_1D(rhs, current_solution, previous_solution, nz, nz_ghost, dz, mpi_info);
    #elif ITERATIVE_SOLVER_TYPE == 1
        gauss_seidel_1D(rhs, current_solution, previous_solution, nz, nz_ghost, dz, mpi_info);
    #endif // ITERATIVE_SOLVER_TYPE
    
    // Updating final_solution
    for (int i = 0; i < nz; i++)
    {
        final_solution[i + nz_ghost] = current_solution[i+1];
    }
    
    deallocate_1D_array(previous_solution);
    deallocate_1D_array(current_solution);
}