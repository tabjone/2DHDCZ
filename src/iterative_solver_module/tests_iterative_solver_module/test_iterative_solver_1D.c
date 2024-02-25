#include "iterative_solver_module/iterative_solver_1D/iterative_solver_1D.h"
#include "array_utilities/array_utilities.h"
#include "MPI_module/MPI_module.h"
#include "global_float_precision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI


#define nz 600

#define epsilon 1.0e-1



void test_iterative_solver_1D()
{
    MPI_Init(NULL, NULL);
    struct MpiInfo *mpi_info = NULL;
    initialize_mpi_info_struct(&mpi_info);

    FLOAT_P *rhs, *analytical_solution, *numerical_solution, *initial_guess;
    allocate_1D_array(&rhs, nz);
    allocate_1D_array(&analytical_solution, nz+4);
    allocate_1D_array(&numerical_solution, nz+4);
    allocate_1D_array(&initial_guess, nz+4);

    FLOAT_P dz = 1.0 / (nz - 1);

    FLOAT_P x;

    for (int i = 0; i < nz; i++)
    {
        x = i * dz;
        rhs[i] = -4.0 * M_PI * M_PI * sin(2.0 * M_PI * x);

    }
    for (int i = 0; i < nz+4; i++)
    {
        x = (i-2) * dz;
        initial_guess[i] = 0.0;
        analytical_solution[i] = sin(2.0 * M_PI * x);
        //printf("x = %f, analytical = %e\n", x, analytical_solution[i]);
    }

    iterative_solver_1D(rhs, numerical_solution, initial_guess, nz, 2, dz, mpi_info);

    for (int i = 0; i < nz+4; i++)
    {        
        printf("%e %e\n", analytical_solution[i], numerical_solution[i]);
    }


    free(mpi_info);
    deallocate_1D_array(rhs);
    deallocate_1D_array(analytical_solution);
    deallocate_1D_array(numerical_solution);
    deallocate_1D_array(initial_guess);

    MPI_Finalize();

}