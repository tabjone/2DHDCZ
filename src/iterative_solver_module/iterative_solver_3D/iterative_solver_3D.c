#include "global_float_precision.h"
#include "global_boundary.h"
#include "global_parameters.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "MPI_module/MPI_module.h"
#include "iterative_solver_3D.h"

void iterative_solver_3D(FLOAT_P ***rhs, FLOAT_P ***final_solution, FLOAT_P ***initial_guess, int nz, int nz_ghost, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, struct MpiInfo *mpi_info)
{
    /*
    Solves the Poisson equation using an iterative method (Jacobi, Gauss-Seidel)

    Parameters
    ----------
    rhs : FLOAT_P ***
        Right hand side of the Poisson equation
    final_solution : FLOAT_P ***
        Array where solution will be stored
    initial_guess : FLOAT_P ***
        Initial guess for the solution
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

    // Initializing previous_solution and current_solution
    FLOAT_P ***current_solution, ***previous_solution;
    allocate_3D_array(&current_solution, nz+2, ny, nx);
    allocate_3D_array(&previous_solution, nz+2, ny, nx);
    
    // Initialize ghost cells
    for (int j = 0; j < ny; j++)
    {
        for (int k = 0; k < nx; k++)
        {
            current_solution[0][j][k] = 0.0;
            current_solution[nz+1][j][k] = 0.0;
        }
    }
    
    // Initialize inside grid
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                current_solution[i+1][j][k] = initial_guess[i+nz_ghost][j][k];
            }
        }
    }

    communicate_3D_ghost_above_below(current_solution, mpi_info, nz, 1, ny, nx);

    // Non-periodic boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (!mpi_info->has_neighbor_below) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    current_solution[1][j][k] = LOWER_PRESSURE_BOUNDARY;
                }
            }
        }
        if (!mpi_info->has_neighbor_above) // Boundary zero
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    current_solution[nz][j][k] = UPPER_PRESSURE_BOUNDARY;
                }
            }
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if ITERATIVE_SOLVER_TYPE == 0
        jacobi_3D(rhs, current_solution, previous_solution, nz, nz_ghost, ny, nx, dz, dy, dx, mpi_info);
    #elif ITERATIVE_SOLVER_TYPE == 1
        gauss_seidel_3D(rhs, current_solution, previous_solution, nz, nz_ghost, ny, nx, dz, dy, dx, mpi_info);
    #endif // ITERATIVE_SOLVER_TYPE
    
    // Updating final_solution
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            final_solution[i + nz_ghost][j] = current_solution[i+1][j];
        }
    }
    
    deallocate_3D_array(previous_solution);
    deallocate_3D_array(current_solution);
}