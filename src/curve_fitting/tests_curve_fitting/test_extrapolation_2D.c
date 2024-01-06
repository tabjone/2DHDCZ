#include "curve_fitting/extrapolation/extrapolation.h"
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "global_float_precision.h"
#include <math.h>
#include <stdio.h>

#define epsilon 1e-3

int test_extrapolate_2D_array_constant_down()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny); 

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }
    
    extrapolate_2D_array_constant_down(array, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] - 1.0) > epsilon ||
            fabs(array[1][j] - 1.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] - 0.0) > epsilon ||
            fabs(array[8][j] - 0.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

int test_extrapolate_2D_array_constant_up()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny); 

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }
    
    extrapolate_2D_array_constant_up(array, nz_full, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] - 0.0) > epsilon ||
            fabs(array[1][j] - 0.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] - 5.0) > epsilon ||
            fabs(array[8][j] - 5.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

int test_extrapolate_2D_array_antisymmetric_down()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny); 

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }
    
    extrapolate_2D_array_antisymmetric_down(array, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] + 3.0) > epsilon ||
            fabs(array[1][j] + 2.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] - 0.0) > epsilon ||
            fabs(array[8][j] - 0.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

int test_extrapolate_2D_array_antisymmetric_up()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny); 

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }
    
    extrapolate_2D_array_antisymmetric_up(array, nz_full, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] - 0.0) > epsilon ||
            fabs(array[1][j] - 0.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] + 4.0) > epsilon ||
            fabs(array[8][j] + 3.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

int test_extrapolate_2D_array_symmetric_down()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny);

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }

    extrapolate_2D_array_symmetric_down(array, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] - 3.0) > epsilon ||
            fabs(array[1][j] - 2.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] - 0.0) > epsilon ||
            fabs(array[8][j] - 0.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

int test_extrapolate_2D_array_symmetric_up()
{
    // Return 0 if test passes, 1 if test fails

    // nz = 5, nz_ghost = 2, nz_full = 9
    int nz_ghost = 2;
    int nz_full = 9;
    int ny = 3;

    FLOAT_P **array;
    allocate_2D_array(&array, nz_full, ny);

    // Initialize ghost cells with 0 and 1->5 in array
    for (int j = 0; j < ny; j++)
    {
        array[0][j] = 0.0;
        array[1][j] = 0.0;
        array[2][j] = 1.0;
        array[3][j] = 2.0;
        array[4][j] = 3.0;
        array[5][j] = 4.0;
        array[6][j] = 5.0;
        array[7][j] = 0.0;
        array[8][j] = 0.0;

    }

    extrapolate_2D_array_symmetric_up(array, nz_full, nz_ghost, ny);

    // Check that ghost cells are extrapolated correctly
    for (int j = 0; j < ny; j++)
    {
        if 
        (
            fabs(array[0][j] - 0.0) > epsilon ||
            fabs(array[1][j] - 0.0) > epsilon ||
            fabs(array[2][j] - 1.0) > epsilon ||
            fabs(array[3][j] - 2.0) > epsilon ||
            fabs(array[4][j] - 3.0) > epsilon ||
            fabs(array[5][j] - 4.0) > epsilon ||
            fabs(array[6][j] - 5.0) > epsilon ||
            fabs(array[7][j] - 4.0) > epsilon ||
            fabs(array[8][j] - 3.0) > epsilon
        )
        {
            deallocate_2D_array(array);
            return 1;
        }
    }
    deallocate_2D_array(array);
    return 0;
}

void test_extrapolation_2D()
{
    if(test_extrapolate_2D_array_constant_down()){
        printf("FAIL: test_extrapolate_2D_array_constant_down\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_constant_down\n");}
    if(test_extrapolate_2D_array_constant_up()){
        printf("FAIL: test_extrapolate_2D_array_constant_up\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_constant_up\n");}
    if(test_extrapolate_2D_array_antisymmetric_down()){
        printf("FAIL: test_extrapolate_2D_array_antisymmetric_down\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_antisymmetric_down\n");}
    if(test_extrapolate_2D_array_antisymmetric_up()){
        printf("FAIL: test_extrapolate_2D_array_antisymmetric_up\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_antisymmetric_up\n");}
    if(test_extrapolate_2D_array_symmetric_down()){
        printf("FAIL: test_extrapolate_2D_array_symmetric_down\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_symmetric_down\n");}
    if(test_extrapolate_2D_array_symmetric_up()){
        printf("FAIL: test_extrapolate_2D_array_symmetric_up\n");}
    else{
        printf("PASS: test_extrapolate_2D_array_symmetric_up\n");}
}