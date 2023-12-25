#include "global_testing.h"
#include <stdio.h>

#include "./spacial_derivatives_module/tests_derivatives_module/tests_derivatives_module.h"
#include "./MPI_module/tests_MPI_module/tests_MPI_module.h"

int main_tests(int argc, char *argv[])
{
    #if TEST_DERIVATIVES_MODULE == 1
        test_derivatives_2D();
    #endif // TESTS_DERIVATIVES_MODULE

    #if TEST_MPI_MODULE == 1
        test_communication();
        //test_mpi_info_struct();
    #endif // TESTS_MPI_MODULE

    return 0;
}