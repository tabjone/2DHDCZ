#include "global_float_precision.h"
#include "global_boundary.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#include "../../MPI_module/MPI_module.h"
#include "iterative_solver_2D.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

void iterative_solver_2D(FLOAT_P **rhs, FLOAT_P **final_solution, FLOAT_P **initial_guess, FLOAT_P, int nz, int nz_ghost, int ny, FLOAT_P dz, FLOAT_P dy, struct MpiInfo *mpi_info)
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
    ny : int
        Number of grid points in y-direction
    dz : FLOAT_P
        Grid spacing in z-direction
    dy : FLOAT_P
        Grid spacing in y-direction
    mpi_info : struct MpiInfo
        Struct containing MPI information
    */

    // Initializing previous_solution and current_solution
    FLOAT_P **current_solution, **previous_solution;
    allocate_2D_array(&current_solution, nz+2, ny);
    allocate_2D_array(&previous_solution, nz+2, ny);
    
    // Initialize ghost cells
    for (int j = 0; j < ny; j++)
    {
        current_solution[0][j] = 0.0;
        current_solution[nz+1][j] = 0.0;
    }
    
    // Initialize inside grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            current_solution[i+1][j] = initial_guess[i+nz_ghost][j];
        }
    }

    communicate_2D_ghost_above_below(current_solution, mpi_info, nz, 1, ny);

    // Non-periodic boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (!mpi_info->has_neighbor_below) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                current_solution[1][j] = LOWER_PRESSURE_BOUNDARY;
            }
        }
        if (!mpi_info->has_neighbor_above) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                current_solution[nz][j] = UPPER_PRESSURE_BOUNDARY;
            }
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if ITERATIVE_SOLVER_TYPE == 0
        jacobi_2D(rhs, current_solution, previous_solution, nz, nz_ghost, ny, dz, dy, mpi_info);
    #elif ITERATIVE_SOLVER_TYPE == 1
        gauss_seidel_2D(rhs, current_solution, previous_solution, nz, nz_ghost, ny, dz, dy, mpi_info);
    #endif // ITERATIVE_SOLVER_TYPE
    
    // Updating final_solution
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            final_solution[i + nz_ghost][j] = current_solution[i+1][j];
        }
    }
    
    deallocate_2D_array(previous_solution);
    deallocate_2D_array(current_solution);
}