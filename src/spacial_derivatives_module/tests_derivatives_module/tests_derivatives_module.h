#ifndef TESTS_DERIVATIVES_MODULE_H
#define TESTS_DERIVATIVES_MODULE_H


int test_central_first_derivative_y_2D();
int test_central_first_derivative_z_2D();
int test_central_second_derivative_y_2D();
int test_central_second_derivative_z_2D();
int test_central_second_derivative_yz_2D();
int test_central_third_derivative_yyz_2D();
int test_central_third_derivative_yzz_2D();
int test_upwind_first_derivative_y_2D();
int test_upwind_first_derivative_z_2D();

int test_derivatives_2D();

#endif // TESTS_DERIVATIVES_MODULE_H