#include "curve_fitting/extrapolation/extrapolation.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "global_float_precision.h"
#include <math.h>
#include <stdio.h>

#define epsilon 1e-3

int test_extrapolate_1D_array_constant_down()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;

    FLOAT_P *array;
    allocate_1D_array(&array, nz_full); 

    // Initialize ghost cells with 0 and 1->5 in array
    array[0] = 0.0;
    array[1] = 0.0;
    array[2] = 1.0;
    array[3] = 2.0;
    array[4] = 3.0;
    array[5] = 4.0;
    array[6] = 5.0;
    array[7] = 0.0;
    array[8] = 0.0;

    extrapolate_1D_array_constant_down(array, nz_ghost);

    // Check that ghost cells are extrapolated correctly
    if 
    (
        fabs(array[0] - 1.0) > epsilon ||
        fabs(array[1] - 1.0) > epsilon ||
        fabs(array[2] - 1.0) > epsilon ||
        fabs(array[3] - 2.0) > epsilon ||
        fabs(array[4] - 3.0) > epsilon ||
        fabs(array[5] - 4.0) > epsilon ||
        fabs(array[6] - 5.0) > epsilon ||
        fabs(array[7] - 0.0) > epsilon ||
        fabs(array[8] - 0.0) > epsilon
    )
    {
        deallocate_1D_array(array);
        return 1;
    }
    deallocate_1D_array(array);
    return 0;
}

int test_extrapolate_1D_array_constant_up()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;

    FLOAT_P *array;
    allocate_1D_array(&array, nz_full); 

    // Initialize ghost cells with 0 and 1->5 in array
    array[0] = 0.0;
    array[1] = 0.0;
    array[2] = 1.0;
    array[3] = 2.0;
    array[4] = 3.0;
    array[5] = 4.0;
    array[6] = 5.0;
    array[7] = 0.0;
    array[8] = 0.0;

    extrapolate_1D_array_constant_up(array, nz_full ,nz_ghost);

    // Check that ghost cells are extrapolated correctly
    if 
    (
        fabs(array[0] - 0.0) > epsilon ||
        fabs(array[1] - 0.0) > epsilon ||
        fabs(array[2] - 1.0) > epsilon ||
        fabs(array[3] - 2.0) > epsilon ||
        fabs(array[4] - 3.0) > epsilon ||
        fabs(array[5] - 4.0) > epsilon ||
        fabs(array[6] - 5.0) > epsilon ||
        fabs(array[7] - 5.0) > epsilon ||
        fabs(array[8] - 5.0) > epsilon
    )
    {
        deallocate_1D_array(array);
        return 1;
    }
    deallocate_1D_array(array);
    return 0;
}

void test_extrapolation_1D()
{
    if(test_extrapolate_1D_array_constant_down()){
        printf("FAIL: test_extrapolate_1D_array_constant_down\n");}
    else{
        printf("PASS: test_extrapolate_1D_array_constant_down\n");}
    if(test_extrapolate_1D_array_constant_up()){
        printf("FAIL: test_extrapolate_1D_array_constant_up\n");}
    else{
        printf("PASS: test_extrapolate_1D_array_constant_up\n");}
}