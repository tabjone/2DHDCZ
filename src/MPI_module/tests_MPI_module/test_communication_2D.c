#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#include "global_float_precision.h"
#include "global_boundary.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "../MPI_module.h"

int test_communication_2D_non_periodic_two_ghosts(struct MpiInfo *mpi_info)
{
    FLOAT_P epsilon = (FLOAT_P)1.0e-9;

    int nz = 5;
    int nz_ghost = 2;
    int ny = 2;

    FLOAT_P **array;
    allocate_2D_array(&array, nz+2*nz_ghost, ny);

    // Ghost cells set to -1.0 and inside the grid we go from 1->5
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = -1.0;
        array[1][j] = -1.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = -1.0;
        array[8][j] = -1.0;
    }

    communicate_2D_ghost_above_below(array, mpi_info, nz, nz_ghost, ny);

    if(mpi_info->rank == 0)
    {   
        // Checking that bottom ghost cells are not changed
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[0][j] + 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
            if (fabs(array[1][j] + 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }
    else
    {
        // Checking that bottom ghost cells are the correct values
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[0][j] - 4.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
            if (fabs(array[1][j] - 5.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }

    if (mpi_info->rank == mpi_info->size - 1)
    {
        // Checking that top ghost cells are not changed
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[7][j] + 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
            if (fabs(array[8][j] + 1.0) > epsilon)
            {
                deallocate_2D_array(array);
            return 1;
            }
        }
    }
    else
    {
        // Checking that top ghost cells are the correct values
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[7][j] - 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
            if (fabs(array[8][j] - 2.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }

    // Checking that inside the grid is unchanged
    for (int j = 0; j < ny; j++)
    {
        if (fabs(array[2][j] - 1.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[3][j] - 2.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[4][j] - 3.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[5][j] - 4.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[6][j] - 5.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
    }

    deallocate_2D_array(array);
    return 0;
}

int test_communication_2D_non_periodic_one_ghost(struct MpiInfo *mpi_info)
{
    FLOAT_P epsilon = (FLOAT_P)1.0e-9;

    int nz = 5;
    int nz_ghost = 1;
    int ny = 2;

    FLOAT_P **array;
    allocate_2D_array(&array, nz+2*nz_ghost, ny);

    // Ghost cells set to -1.0 and inside the grid we go from 1->5
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = -1.0;
        array[1][j] = 1.0;
        array[2][j] = 2.0;
        array[3][j] = 3.0;
        array[4][j] = 4.0;
        array[5][j] = 5.0;
        array[6][j] = -1.0;
    }

    communicate_2D_ghost_above_below(array, mpi_info, nz, nz_ghost, ny);

    if(mpi_info->rank == 0)
    {   
        // Checking that bottom ghost cells are not changed
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[0][j] + 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }
    else
    {
        // Checking that bottom ghost cells are the correct values
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[0][j] - 5.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }

    if (mpi_info->rank == mpi_info->size - 1)
    {
        // Checking that top ghost cells are not changed
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[6][j] + 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }
    else
    {
        // Checking that top ghost cells are the correct values
        for (int j = 0; j < ny; j++)
        {
            if (fabs(array[6][j] - 1.0) > epsilon)
            {
            deallocate_2D_array(array);
            return 1;
            }
        }
    }

    // Checking that inside the grid is unchanged
    for (int j = 0; j < ny; j++)
    {
        if (fabs(array[1][j] - 1.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[2][j] - 2.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[3][j] - 3.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[4][j] - 4.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
        if (fabs(array[5][j] - 5.0) > epsilon)
        {
            deallocate_2D_array(array);
            return 1;
        }
    }

    deallocate_2D_array(array);
    return 0;
}

void test_communication_2D()
{
    
    MPI_Init(NULL, NULL);
    struct MpiInfo *mpi_info = NULL;
    initialize_mpi_info_struct(&mpi_info);

    int my_test, test;

    #if PERIODIC_BOUNDARY_VERTICAL == 0
        my_test = test_communication_2D_non_periodic_one_ghost(mpi_info)
        + test_communication_2D_non_periodic_two_ghosts(mpi_info);
    
    #else
        //my_test = ....
    #endif // PERIODIC_BOUNDARY_VERTICAL

    MPI_Allreduce(&my_test, &test, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_info->rank == 0)
    {
        if (test!=0)
        {printf("FAIL: 2D communication.\n");}
        else
        {printf("PASS: 2D communication.\n");}
    }


    free(mpi_info);
    MPI_Finalize();
}