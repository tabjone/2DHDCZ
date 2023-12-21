#include <stdlib.h>
#include "global_float_precision.h"

void allocate_3D_array(FLOAT_P ****array_ptr, int nz, int ny, int nx) 
{
    /*
    Allocates memory for a 3D contiguous array.

    Parameters
    ----------
    array_ptr : FLOAT_P****
        Pointer to the array to allocate.
    nz : int
        Number of layers in the 3D array (depth).
    ny : int
        Number of rows in each 2D slice of the 3D array.
    nx : int
        Number of columns in each 2D slice of the 3D array.
    */

    // Allocate the primary array containing pointers to each 2D slice.
    *array_ptr = malloc(nz * sizeof(***array_ptr));
    if (*array_ptr == NULL) {
        // handle malloc failure
        return;
    }

    // Allocate the array containing pointers to each row of each 2D slice.
    (*array_ptr)[0] = malloc(nz*ny * sizeof(**array_ptr[0]));
    if ((*array_ptr)[0] == NULL) {
        // handle malloc failure
        free(*array_ptr);
        return;
    }

    // Allocate the contiguous block of memory for the actual 3D array.
    (*array_ptr)[0][0] = malloc(nz*ny*nx * sizeof(*array_ptr[0][0])); // Contiguous
    if ((*array_ptr)[0][0] == NULL) {
        // handle malloc failure
        free((*array_ptr)[0]);
        free(*array_ptr);
        return;
    }

    // Set up the pointers into the contiguous memory.
    for (int i = 1; i < nz; i++) {
        (*array_ptr)[i] = &(*array_ptr)[0][i*ny];
    }

    for (int j = 1; j < nz*ny; j++) {
        (*array_ptr)[0][j] = &(*array_ptr)[0][0][j*nx];
    }
}
