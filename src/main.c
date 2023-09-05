#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

#include "functions.h"
#include "global_parameters.h"

int main(int argc, char *argv[])
{
    #if MHD == 0
        #if DIMENSIONS == 1
            printf("Not implemented yet\n");
        #elif DIMENSIONS == 2
            main_hd_2D(argc, argv);
        #elif DIMENSIONS == 3
            printf("Not implemented yet\n");
        #endif
    #elif MHD == 1
        #if DIMENSIONS == 1
            printf("Not implemented yet\n");
        #elif DIMENSIONS == 2
            printf("Not implemented yet\n");
        #elif DIMENSIONS == 3
            printf("Not implemented yet\n");
        #endif
    #endif
}
