#include <stdlib.h>

void deallocate_1D_array(double *array_ptr) 
{
    /*
    Deallocates memory for a 1D array.

    Parameters
    ----------
    array_ptr : double**
        Pointer to the 1D array to deallocate.
    */

    free(array_ptr);  // Free the actual 1D array
    array_ptr = NULL; // Set the pointer to NULL to avoid dangling pointer
}
