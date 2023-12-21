#include <stdio.h>
#include "math.h"

#include "../derivatives_2D/derivatives_2D.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI


#define Lz (FLOAT_P)1.0e8
#define Ly (FLOAT_P)1.0e8
#define nz 100
#define ny 100
#define epsilon (FLOAT_P)1.0e-2

int test_central_first_derivative_y_2D() 
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P one_over_2dy = 1.0 / (2.0 * dy);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = 2.0*M_PI/Ly*sin(4.0*M_PI*z/Lz)*cos(2*M_PI*y/Ly);
            computed_value = central_first_derivative_y_2D(f, i, j, ny, one_over_2dy);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_first_derivative_z_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = 4.0*M_PI/Lz*cos(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Lz);
            computed_value = central_first_derivative_z_2D(f, i, j, one_over_2dz);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                printf("i: %d, j: %d\n", i, j);
                printf("expected_value: %f, computed_value: %f\n", expected_value, computed_value);
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_second_derivative_y_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_dydy = 1.0 / (dy * dy);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = -4.0*M_PI*M_PI/Ly/Ly*sin(4.0*M_PI*z/Lz)*sin(2.0*M_PI*y/Ly);
            computed_value = central_second_derivative_y_2D(f, i, j, ny, one_over_dydy);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_second_derivative_z_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_dzdz = 1.0 / (dz * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = -16.0*M_PI*M_PI/Lz/Lz*sin(4.0*M_PI*z/Lz)*sin(2.0*M_PI*y/Ly);
            computed_value = central_second_derivative_z_2D(f, i, j, one_over_dzdz);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_second_derivative_yz_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_4dydz = 1.0 / (4.0 * dy * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = 8.0*M_PI*M_PI/Lz/Ly*cos(4.0*M_PI*z/Lz)*cos(2.0*M_PI*y/Ly);
            computed_value = central_second_derivative_yz_2D(f, i, j, ny, one_over_4dydz);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_third_derivative_yyz_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_8dydydz = 1.0 / (8.0 * dy * dy * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = -16.0*M_PI*M_PI*M_PI/Ly/Ly/Lz*cos(4.0*M_PI*z/Lz)*sin(2.0*M_PI*y/Ly);
            computed_value = central_third_derivative_yyz_2D(f, i, j, ny, one_over_8dydydz);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_central_third_derivative_yzz_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;

    allocate_2D_array(&f, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_8dydzdz = 1.0 / (8.0 * dy * dz * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 2; i < nz-2; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            expected_value = -32.0*M_PI*M_PI*M_PI/Lz/Lz/Ly*sin(4.0*M_PI*z/Lz)*cos(2.0*M_PI*y/Ly);
            computed_value = central_third_derivative_yzz_2D(f, i, j, ny, one_over_8dydzdz);
            if (fabs(expected_value - computed_value) > epsilon)
            {
                deallocate_2D_array(f);
                return 1;
            }            
        }
    }

    deallocate_2D_array(f);
    return 0;
}

int test_upwind_first_derivative_y_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;
    FLOAT_P **velocity_positive, **velocity_negative;

    allocate_2D_array(&f, nz, ny);
    allocate_2D_array(&velocity_positive, nz, ny);
    allocate_2D_array(&velocity_negative, nz, ny);

    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P one_over_dy = 1.0 / dy;
    FLOAT_P one_over_2dy = 1.0 / (2.0 * dy);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 1; j < ny-1; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
            velocity_positive[i][j] = 1.0;
            velocity_negative[i][j] = -1.0;
        }
    }

    FLOAT_P expected_value, computed_value_positive, computed_value_negative;

    // Calculating derivative
    for (int i = 0; i < nz; i++) {
        for (int j = 0; j < ny; j++) {
            expected_value = 2.0*M_PI/Ly*sin(4.0*M_PI*z/Lz)*cos(2*M_PI*y/Ly);
            computed_value_positive = upwind_first_derivative_y_2D(f, velocity_positive, i, j, ny, one_over_dy, one_over_2dy);
            computed_value_negative = upwind_first_derivative_y_2D(f, velocity_negative, i, j, ny, one_over_dy, one_over_2dy);
            if (fabs(expected_value - computed_value_positive) > epsilon || fabs(expected_value - computed_value_negative) > epsilon)
            {
                deallocate_2D_array(f);
                deallocate_2D_array(velocity_positive);
                deallocate_2D_array(velocity_negative);
                return 1;
            }
        }
    }

    deallocate_2D_array(f);
    deallocate_2D_array(velocity_positive);
    deallocate_2D_array(velocity_negative);
    return 0;
}

int test_upwind_first_derivative_z_2D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P **f;
    FLOAT_P **velocity_positive, **velocity_negative;

    allocate_2D_array(&f, nz, ny);
    allocate_2D_array(&velocity_positive, nz, ny);
    allocate_2D_array(&velocity_negative, nz, ny);

    FLOAT_P dy = 1.0 / (ny-1);
    FLOAT_P dz = 1.0 / (nz-1);
    FLOAT_P one_over_dz = 1.0 / dz;
    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P z, y;
    for (int i = 0; i < nz; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            y = j * dy;
            f[i][j] = sin(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
            velocity_positive[i][j] = 1.0;
            velocity_negative[i][j] = -1.0;
        }
    }

    FLOAT_P expected_value, computed_value_positive, computed_value_negative;

    // Calculating derivative
    for (int i = 2; i < nz-2; i++) {
        z = i * dz;
        for (int j = 0; j < ny; j++) {
            expected_value = 4.0*M_PI/Lz*cos(4.0*M_PI*z/Lz)*sin(2*M_PI*y/Ly);
            computed_value_positive = upwind_first_derivative_z_2D(f, velocity_positive, i, j, one_over_dz, one_over_2dz);
            computed_value_negative = upwind_first_derivative_z_2D(f, velocity_negative, i, j, one_over_dz, one_over_2dz);
            if (fabs(expected_value - computed_value_positive) > epsilon || fabs(expected_value - computed_value_negative) > epsilon)
            {
                deallocate_2D_array(f);
                deallocate_2D_array(velocity_positive);
                deallocate_2D_array(velocity_negative);
                return 1;
            }
        }
    }

    deallocate_2D_array(f);
    deallocate_2D_array(velocity_positive);
    deallocate_2D_array(velocity_negative);
    return 0;
}


int test_derivatives_2D()
{
    if(test_central_first_derivative_y_2D()){
        printf("FAIL: Central first derivative in y-direction.\n");}
    else{
        printf("PASS: Central first derivative in y-direction.\n");}
    if(test_central_first_derivative_z_2D()){
        printf("FAIL: Central first derivative in z-direction.\n");}
    else{
        printf("PASS: Central first derivative in z-direction.\n");}
    
    if(test_central_second_derivative_y_2D()){
        printf("FAIL: Central second derivative in y-direction.\n");}
    else{
        printf("PASS: Central second derivative in y-direction.\n");}
    if(test_central_second_derivative_z_2D()){
        printf("FAIL: Central second derivative in z-direction.\n");}
    else{
        printf("PASS: Central second derivative in z-direction.\n");}
    if(test_central_second_derivative_yz_2D()){
        printf("FAIL: Central second derivative in yz-direction.\n");}
    else{
        printf("PASS: Central second derivative in yz-direction.\n");}
    if(test_central_third_derivative_yyz_2D()){
        printf("FAIL: Central third derivative in yyz-direction.\n");}
    else{
        printf("PASS: Central third derivative in yyz-direction.\n");}
    if(test_central_third_derivative_yzz_2D()){
        printf("FAIL: Central third derivative in yzz-direction.\n");}
    else{
        printf("PASS: Central third derivative in yzz-direction.\n");}
    if(test_upwind_first_derivative_y_2D()){
        printf("FAIL: Upwind first derivative in y-direction.\n");}
    else{
        printf("PASS: Upwind first derivative in y-direction.\n");}
    if(test_upwind_first_derivative_z_2D()){
        printf("FAIL: Upwind first derivative in z-direction.\n");}
    else{
        printf("PASS: Upwind first derivative in z-direction.\n");}
    return 0;
}