#include <stdio.h>
#include "math.h"
#include "global_float_precision.h"

#include "../derivatives_1D/derivatives_1D.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define Lz (FLOAT_P)1.0e1
#define nz 200
#define epsilon (FLOAT_P)1.0e-2

int test_central_first_derivative_z_1D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P *f;

    allocate_1D_array(&f, nz);

    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        f[i] = sin(4.0*M_PI*z/Lz);
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        expected_value = 4.0*M_PI/Lz*cos(4.0*M_PI*z/Lz);
        computed_value = central_first_derivative_z_1D(f, i, one_over_2dz);
        if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
        {
            deallocate_1D_array(f);
            return 1;
        }            
    }

    deallocate_1D_array(f);
    return 0;
}

int test_central_second_derivative_z_1D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P *f;

    allocate_1D_array(&f, nz);

    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_dzdz = 1.0 / (dz * dz);

    FLOAT_P z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        f[i] = sin(4.0*M_PI*z/Lz);
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        expected_value = -16.0*M_PI*M_PI/(Lz*Lz)*sin(4.0*M_PI*z/Lz);
        computed_value = central_second_derivative_z_1D(f, i, one_over_dzdz);
        if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
        {
            deallocate_1D_array(f);
            return 1;
        }            
    }
    deallocate_1D_array(f);
    return 0;
}

int test_upwind_first_derivative_z_1D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P *f;
    FLOAT_P *velocity_positive, *velocity_negative;

    allocate_1D_array(&f, nz);
    allocate_1D_array(&velocity_positive, nz);
    allocate_1D_array(&velocity_negative, nz);

    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_dz = 1.0 / dz;
    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        f[i] = sin(4.0*M_PI*z/Lz);
        velocity_negative[i] = -1.0;
        velocity_positive[i] = 1.0;
    }

    FLOAT_P expected_value, computed_value_positive, computed_value_negative;

    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {
        z = i * dz;
        expected_value = 4.0*M_PI/Lz*cos(4.0*M_PI*z/Lz);
        computed_value_negative = upwind_first_derivative_z_1D(f, velocity_negative, i, one_over_dz, one_over_2dz);
        computed_value_positive = upwind_first_derivative_z_1D(f, velocity_positive, i, one_over_dz, one_over_2dz);

        if (fabs(expected_value - computed_value_negative)/fabs(expected_value) > epsilon || fabs(expected_value - computed_value_positive)/fabs(expected_value) > epsilon)
        {
            deallocate_1D_array(f);
            deallocate_1D_array(velocity_positive);
            deallocate_1D_array(velocity_negative);
            return 1;
        }            
    }

    deallocate_1D_array(f);
    deallocate_1D_array(velocity_positive);
    deallocate_1D_array(velocity_negative);
    return 0;
}

int test_derivatives_1D()
{
    if(test_central_first_derivative_z_1D()){
        printf("FAIL: test_central_first_derivative_z_1D.\n");}
    else{
        printf("PASS: test_central_first_derivative_z_1D.\n");}
    if(test_central_second_derivative_z_1D()){
        printf("FAIL: test_central_second_derivative_z_1D.\n");}
    else{
        printf("PASS: test_central_second_derivative_z_1D.\n");}
    if(test_upwind_first_derivative_z_1D()){
        printf("FAIL: test_upwind_first_derivative_z_1D.\n");}
    else{
        printf("PASS: test_upwind_first_derivative_z_1D.\n");}
    return 0;
}