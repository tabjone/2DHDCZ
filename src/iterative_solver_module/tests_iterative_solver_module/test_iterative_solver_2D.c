#include "iterative_solver_module/iterative_solver_2D/iterative_solver_2D.h"
#include "array_utilities/array_utilities.h"
#include "MPI_module/MPI_module.h"
#include "global_float_precision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI


#define nz 10
#define ny 10

#define epsilon 1.0e-1

// Testing for analytical solution p = sin(2*pi*x)*cos(2*pi*y)
// d^2p/dx^2 + d^2p/dy^2 = - 8 pi^2 sin(2*pi*x)*cos(2*pi*y)

void test_iterative_solver_2D()
{
    MPI_Init(NULL, NULL);
    struct MpiInfo *mpi_info = NULL;
    initialize_mpi_info_struct(&mpi_info);

    FLOAT_P **rhs, **analytical_solution, **numerical_solution, **initial_guess;
    allocate_2D_array(&rhs, nz, ny);
    allocate_2D_array(&analytical_solution, nz+4, ny);
    allocate_2D_array(&numerical_solution, nz+4, ny);
    allocate_2D_array(&initial_guess, nz+4, ny);

    FLOAT_P dz = 1.0 / (nz - 1);
    FLOAT_P dy = 1.0 / (ny - 1);

    FLOAT_P x, y;

    for (int i = 0; i < nz; i++)
    {
        x = i * dz;
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            rhs[i][j] = -8.0 * M_PI * M_PI * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * y);
        }
    }
    for (int i = 0; i < nz+4; i++)
    {
        x = (i-2) * dz;
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            initial_guess[i][j] = 0.0;
            analytical_solution[i][j] = sin(2.0 * M_PI * x) * cos(2.0 * M_PI * y);
        }    
    }

    iterative_solver_2D(rhs, numerical_solution, initial_guess, nz, 2, ny, dz, dy, mpi_info);

    FLOAT_P abs_difference;

    printf("[");
    for (int i = 0; i < nz+4; i++)
    {
        printf("[");        
        for (int j = 0; j < ny; j++)
        {
            abs_difference = fabs(analytical_solution[i][j] - numerical_solution[i][j]);
            printf("%e, ", numerical_solution[i][j]);
        }
        printf("], \n");
    }
    printf("]\n");


    free(mpi_info);
    deallocate_2D_array(rhs);
    deallocate_2D_array(analytical_solution);
    deallocate_2D_array(numerical_solution);
    deallocate_2D_array(initial_guess);

    MPI_Finalize();

}