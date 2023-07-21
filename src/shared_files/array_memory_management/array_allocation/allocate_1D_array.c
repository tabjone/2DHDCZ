#include <stdlib.h>

void allocate_1D_array(double **array_ptr, int nx) 
{
    /*
    Allocates memory for a 1D array.

    Parameters
    ----------
    array_ptr : double**
        Pointer to the array to allocate.
    nx : int
        Size of the array.
    */
   
    // Allocate the primary array.
    *array_ptr = malloc(nx * sizeof(**array_ptr));
    if (*array_ptr == NULL) {
        return;  // Failed to allocate memory
    }
}