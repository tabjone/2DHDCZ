#include <stdlib.h>
#include "global_parameters.h"


void deallocate_2D_array(FLOAT_P **array_ptr) 
{
    /*
    Deallocates memory for a 2D contiguous array.

    Parameters
    ----------
    array_ptr : FLOAT_P**
        Pointer to the array to deallocate.
    */

    free(array_ptr[0]);  // Free the actual 2D array
    free(array_ptr);     // Free the primary array
    array_ptr = NULL;      // Set the pointer to NULL to avoid dangling pointer
}