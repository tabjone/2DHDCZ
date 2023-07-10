#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"


#include "functions.h"

int main(int argc, char *argv[])
{   
    double ***threeD_test;
    allocate_3D_array(&threeD_test, 10, 10, 10);
    deallocate_3D_array(threeD_test);

    double **twoD_test;
    allocate_2D_array(&twoD_test, 10, 10);
    deallocate_2D_array(twoD_test);

    double *oneD_test;
    allocate_1D_array(&oneD_test, 10);
    deallocate_1D_array(oneD_test);
}
