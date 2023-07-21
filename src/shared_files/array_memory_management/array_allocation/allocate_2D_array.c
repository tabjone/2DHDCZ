#include <stdlib.h>

void allocate_2D_array(double ***array_ptr, int ny, int nx) 
{
    /*
    Allocates memory for a 2D contiguous array.

    Parameters
    ----------
    array_ptr : double***
        Pointer to the array to allocate.
    ny : int
        Number of rows in the array.
    nx : int
        Number of columns in the array.
    */
   
    // Allocate the primary array.
    *array_ptr = malloc(ny * sizeof(**array_ptr));
    if (*array_ptr == NULL) {
        return;  // Failed to allocate memory
    }
    
    // Allocate the actual 2D array.
    (*array_ptr)[0] = malloc(ny * nx * sizeof *(*array_ptr[0]));
    if ((*array_ptr)[0] == NULL) {
        free(*array_ptr);  // Failed to allocate memory, so clean up.
        return;
    }

    // Set up the pointers into the contiguous memory.
    for (int i = 1; i < ny; i++) {
        (*array_ptr)[i] = &((*array_ptr)[0][i * nx]);
    }
}