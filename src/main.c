#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "functions.h"

int main(int argc, char *argv[])
{
    char dir_name[100];
    struct stat st = {0}; // Initialize the stat structure

    // Use snprintf to create the directory path with the base directory and RUN_NAME
    snprintf(dir_name, sizeof(dir_name), "data/%s", RUN_NAME);

    // Check if directory exists
    if (stat(dir_name, &st) == -1) {
        // If directory doesn't exist, create it with permissions set to rwxr-xr-x
        if (mkdir(dir_name, 0755) == -1) {
            perror("Error creating directory");
            return 1;
        }
    }

    #if DIMENSIONS == 1
        printf("Not implemented yet\n");
    #elif DIMENSIONS == 2
        main_2D(argc, argv);
    #elif DIMENSIONS == 3
        printf("Not implemented yet\n");
    #endif // DIMENSIONS
}
