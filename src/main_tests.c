#include "global_testing.h"

#include "./spacial_derivatives_module/tests_derivatives_module/tests_derivatives_module.h"
#include "./MPI_module/tests_MPI_module/tests_MPI_module.h"
#include "./curve_fitting/tests_curve_fitting/tests_curve_fitting.h"
#include "./iterative_solver_module/tests_iterative_solver_module/tests_iterative_solver_module.h"
#include "solver/rhs_functions/tests_rhs_functions/tests_rhs_functions.h"


int main_tests(int argc, char *argv[])
{
    #if TEST_DERIVATIVES_MODULE == 1
        test_derivatives_1D();
        test_derivatives_2D();
        test_derivatives_3D();
    #endif // TESTS_DERIVATIVES_MODULE

    #if TEST_MPI_MODULE == 1
        //test_communication_1D();
        test_communication_2D();
        //test_mpi_info_struct();
    #endif // TESTS_MPI_MODULE

    #if TEST_CURVE_FITTING == 1
        test_extrapolation_1D();
        test_extrapolation_2D();
        test_extrapolation_3D();
    #endif // TESTS_CURVE_FITTING

    #if TEST_ITERATIVE_SOLVER == 1
        //test_iterative_solver_1D();
        test_iterative_solver_2D();
    #endif // TESTS_ITERATIVE_SOLVER

    #if TEST_RHS_FUNCTIONS == 1
        test_rhs_functions_1D();
    #endif // TESTS_RHS_FUNCTIONS

    return 0;
}