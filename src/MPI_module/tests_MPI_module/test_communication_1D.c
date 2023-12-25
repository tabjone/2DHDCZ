#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#include "global_float_precision.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#include "../mpi_info_struct.h"
#include "../communication/communication.h"

int test_communication_1D(struct MpiInfo *mpi_info, int nz, int nz_ghost, FLOAT_P epsilon)
{
    FLOAT_P *array;
    allocate_1D_array(&array, nz+2*nz_ghost);

    // Ghost cells set to -1.0 and inside the grid we go from 0->nz-1
    for (int i = 0; i < nz_ghost; i++)
    {
        array[i] = -1.0;
        array[nz + nz_ghost + i] = -1.0;
    }
    for (int i = nz_ghost; i < nz + nz_ghost; i++)
    {
        array[i] = (FLOAT_P)(i - nz_ghost);
    }

    communicate_1D_ghost_above_below(array, mpi_info, nz, nz_ghost);

    if(mpi_info->rank ==0)
    {
        // Checking that bottom ghost cells are not changed
        for (int i = 0; i < nz_ghost; i++)
        if (fabs(array[i] + 1) > epsilon)
        {
            deallocate_1D_array(array);
            return 1;
        }
    }
    else
    {
        // Checking that bottom ghost cells hold the correct values
        for (int i = 0; i < nz_ghost; i++)
        if (fabs(array[i] - (FLOAT_P)(i + nz)) > epsilon)
        {
            deallocate_1D_array(array);
            return 1;
    }
        // Checking that inside of grid and boundaries are not changed
        for (int i = nz_ghost; i < nz + nz_ghost; i++)
        if (fabs(array[i] - (FLOAT_P)(i - nz_ghost)) > epsilon)
        {
            deallocate_1D_array(array);
            return 1;
        }
    }
    if (mpi_info->rank == mpi_info->size - 1)
    {
        // Checking that top ghost cells are not changed
        for (int i = nz + nz_ghost; i < nz + 2*nz_ghost; i++)
        if (fabs(array[i] + 1) > epsilon)
        {
            deallocate_1D_array(array);
            return 1;
        }
    }
    else
    {
        // Checking that top ghost cells hold the correct values
        for (int i = nz + nz_ghost; i < nz + 2*nz_ghost; i++)
        if (fabs(array[i] - (FLOAT_P)(i - nz - nz_ghost)) > epsilon)
        {
            deallocate_1D_array(array);
            return 1;
        }
    }

    deallocate_1D_array(array);
    return 0;
}

void test_communication_1D()
{
    int nz = 6;
    int nz_ghost = 3;
    FLOAT_P epsilon = (FLOAT_P)1.0e-9;

    MPI_Init(NULL, NULL);
    struct MpiInfo *mpi_info = NULL;
    initialize_mpi_info_struct(&mpi_info);

    int my_test, test;

    #if PERIODIC_BOUNDARY_VERTICAL == 0
        my_test = test_communication_1D(mpi_info, nz, nz_ghost, epsilon);
    
    #else
        my_test = ....
    #endif // PERIODIC_BOUNDARY_VERTICAL

    MPI_Allreduce(&my_test, &test, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (mpi_info->rank == 0)
    {
        if (test)
        {printf("FAIL: 1D communication.\n");}
        else
        {printf("PASS: 1D communication.\n");}
    }


    free(mpi_info);
    MPI_Finalize();
}