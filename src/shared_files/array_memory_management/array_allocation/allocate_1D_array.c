#include <stdlib.h>
#include "global_parameters.h"


void allocate_1D_array(FLOAT_P **array_ptr, int nx) 
{
    /*
    Allocates memory for a 1D array.

    Parameters
    ----------
    array_ptr : FLOAT_P**
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