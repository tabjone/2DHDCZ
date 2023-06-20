#include <stdlib.h>

void deallocate_2D_array(double **array_ptr) {
    free(array_ptr[0]);  // Free the actual 2D array
    free(array_ptr);     // Free the primary array
}